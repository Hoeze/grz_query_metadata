from __future__ import annotations

import collections
import json

import pytest

from grz_query_metadata import survey as mod


def counters() -> mod.Counters:
    return collections.defaultdict(collections.Counter)


class TestWalk:
    def test_plain_path(self):
        assert mod.collect({"a": {"b": 1}}, "a/b") == [1]

    def test_descends_into_arrays(self):
        meta = {"donors": [{"labData": [{"x": 1}, {"x": 2}]}, {"labData": [{"x": 3}]}]}
        assert mod.collect(meta, "donors[]/labData[]/x") == [1, 2, 3]

    def test_missing_key_yields_nothing(self):
        assert mod.collect({"a": {}}, "a/b/c") == []

    def test_null_values_are_dropped(self):
        assert mod.collect({"donors": [{"gender": None}, {"gender": "male"}]}, "donors[]/gender") == ["male"]

    def test_array_segment_over_non_list(self):
        assert mod.collect({"donors": {"gender": "male"}}, "donors[]/gender") == []


class TestCountSubmission:
    def test_counts_enum_and_freetext(self, submission):
        c = counters()
        mod.count_submission(submission(), c)
        assert c["labData.libraryType"]["wgs"] == 1
        assert c["files.fileType"]["fastq"] == 2
        assert c["labData.sequencerModel"]["Illumina NovaSeq 6000"] == 1

    @pytest.mark.parametrize(
        "tissue_type_id,expected",
        [
            ("BTO:0000089", "yes"),
            ("BTO:89", "no"),
            ("whole blood", "no"),
            ("BTO:0000089\n", "no"),  # $ would accept the trailing newline; \Z must not
        ],
    )
    def test_derived_bto_format(self, tissue_type_id, expected, submission):
        c = counters()
        mod.count_submission(submission(tissue_type_id=tissue_type_id), c)
        assert c["tissueTypeId_is_BTO_format"][expected] == 1

    @pytest.mark.parametrize(
        "meta",
        [
            {"donors": [None]},
            {"donors": "not a list"},
            {"donors": [{"labData": {"not": "a list"}}]},
            {"donors": [{"labData": [None]}]},
        ],
    )
    def test_malformed_shapes_are_skipped_not_fatal(self, meta):
        # One odd row must never abort a whole production run.
        c = counters()
        mod.count_submission(meta, c)
        assert sum(c["tissueTypeId_is_BTO_format"].values()) == 0


class TestInitialSubmissions:
    """Duplicate initial submissions per Leistungserbringer."""

    def fold(self, *metas, qc_passed: bool | None = True) -> mod.Counters:
        """Fold submissions in, taking each case id from the metadata the way the
        `pseudonym` column would hold it on an unredacted row."""
        c = counters()
        initials = mod.InitialSubmissions()
        for meta in metas:
            initials.add(meta, qc_passed, (meta.get("submission") or {}).get("localCaseId"))
        initials.write(c)
        return c

    def test_one_initial_per_case_is_no_duplicate(self, submission):
        c = self.fold(submission(local_case_id="a"), submission(local_case_id="b"))
        # One LE, nothing duplicated — it still has to show up in the distribution.
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter({"0": 1})

    def test_a_second_initial_for_the_same_case_is_one_duplicate(self, submission):
        c = self.fold(submission(local_case_id="a"), submission(local_case_id="a"))
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter({"1": 1})

    @pytest.mark.parametrize("sent,duplicates", [(3, "2"), (7, "6"), (30, "29")])
    def test_the_row_label_is_the_exact_count(self, sent, duplicates, submission):
        c = self.fold(*[submission(local_case_id="a")] * sent)
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter({duplicates: 1})

    def test_duplicates_of_two_cases_add_up_for_the_same_LE(self, submission):
        c = self.fold(
            *[submission(local_case_id="a")] * 2,
            *[submission(local_case_id="b")] * 3,
        )
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter({"3": 1})

    def test_the_same_case_id_at_two_LEs_is_not_a_duplicate(self, submission):
        # localCaseId is the LE's own numbering, so it only means anything
        # within one LE; two of them may well pick "case-1".
        c = self.fold(
            submission(submitter_id="260000001", local_case_id="a"),
            submission(submitter_id="260000002", local_case_id="a"),
        )
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter({"0": 2})

    def test_each_LE_is_counted_separately(self, submission):
        c = self.fold(
            *[submission(submitter_id="260000001", local_case_id="a")] * 2,
            submission(submitter_id="260000002", local_case_id="a"),
        )
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter({"1": 1, "0": 1})

    @pytest.mark.parametrize("submission_type", ["followup", "addition", "correction"])
    def test_only_initial_submissions_count(self, submission_type, submission):
        # Sending the same case again is what these types are FOR.
        c = self.fold(*[submission(submission_type=submission_type, local_case_id="a")] * 2)
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter()
        assert c[mod.INITIALS_CHECKED] == collections.Counter(
            {"yes": 0, mod.NOT_QC_PASSED: 0, mod.REDACTED_CASE_ID: 0, mod.NOT_ATTRIBUTABLE: 0}
        )

    @pytest.mark.parametrize("missing", [{"submitter_id": None}, {"local_case_id": None}])
    def test_an_initial_missing_an_identifier_is_counted_apart(self, missing, submission):
        # It cannot be grouped, so it must not quietly lower the duplicate count.
        c = self.fold(submission(local_case_id="a"), submission(**missing))
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter({"0": 1})
        assert c[mod.INITIALS_CHECKED] == collections.Counter(
            {"yes": 1, mod.NOT_QC_PASSED: 0, mod.REDACTED_CASE_ID: 0, mod.NOT_ATTRIBUTABLE: 1}
        )

    @pytest.mark.parametrize("redacted", ["REDACTED_LOCAL_CASE_ID", ""])
    def test_a_redacted_case_id_groups_nothing(self, redacted, submission):
        # grz-db redacts localCaseId inside the stored metadata document. Keying
        # on the placeholder would report every redacted submission of an LE as
        # a duplicate of every other.
        initials = mod.InitialSubmissions()
        for _ in range(4):
            initials.add(submission(local_case_id=redacted), True, redacted)
        c = counters()
        initials.write(c)
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter()
        assert c[mod.INITIALS_CHECKED][mod.REDACTED_CASE_ID] == 4

    def test_the_column_wins_over_the_redacted_metadata(self, submission):
        # The document says REDACTED_LOCAL_CASE_ID for both; the column knows
        # they are two different cases, so neither is a duplicate.
        initials = mod.InitialSubmissions()
        for case in ("real-1", "real-2"):
            initials.add(submission(local_case_id="REDACTED_LOCAL_CASE_ID"), True, case)
        c = counters()
        initials.write(c)
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter({"0": 1})
        assert c[mod.INITIALS_CHECKED]["yes"] == 2

    @pytest.mark.parametrize("meta", [{}, {"submission": None}, {"submission": "not a dict"}])
    def test_malformed_shapes_are_skipped_not_fatal(self, meta):
        initials = mod.InitialSubmissions()
        initials.add(meta, True, "a")
        c = counters()
        initials.write(c)
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter()

    @pytest.mark.parametrize("qc_passed", [False, None])
    def test_a_submission_that_did_not_pass_qc_is_not_counted(self, qc_passed, submission):
        # A rejected submission was MEANT to be sent again; counting the
        # replacement would report the process working as a fault. None is
        # "not decided yet", which is not a pass either.
        c = self.fold(*[submission(local_case_id="a")] * 2, qc_passed=qc_passed)
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter()
        assert c[mod.INITIALS_CHECKED][mod.NOT_QC_PASSED] == 2

    def test_only_the_qc_passed_ones_of_a_case_are_compared(self, submission):
        initials = mod.InitialSubmissions()
        initials.add(submission(local_case_id="a"), False, "a")  # rejected
        initials.add(submission(local_case_id="a"), True, "a")  # the accepted retry
        c = counters()
        initials.write(c)
        # One accepted initial submission for the case: nothing was duplicated.
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter({"0": 1})
        assert c[mod.INITIALS_CHECKED] == collections.Counter(
            {"yes": 1, mod.NOT_QC_PASSED: 1, mod.REDACTED_CASE_ID: 0, mod.NOT_ATTRIBUTABLE: 0}
        )

    def test_two_qc_passed_initials_for_one_case_are_a_duplicate(self, submission):
        initials = mod.InitialSubmissions()
        for qc in (True, False, True):
            initials.add(submission(local_case_id="a"), qc, "a")
        c = counters()
        initials.write(c)
        assert c[mod.DUPLICATE_INITIALS] == collections.Counter({"1": 1})

    def test_duplicated_cases_names_them_worst_first(self, submission):
        initials = mod.InitialSubmissions()
        for meta in (
            *[submission(submitter_id="LE1", local_case_id="a")] * 3,
            *[submission(submitter_id="LE1", local_case_id="b")] * 2,
            submission(submitter_id="LE1", local_case_id="c"),  # sent once: not duplicated
            *[submission(submitter_id="LE2", local_case_id="a")] * 2,
        ):
            initials.add(meta, True, (meta["submission"] or {}).get("localCaseId"))
        assert initials.duplicated_cases() == {
            "LE1": [("a", 3), ("b", 2)],
            "LE2": [("a", 2)],
        }

    def test_no_identifier_ever_reaches_the_counters(self, submission):
        c = self.fold(submission(submitter_id="260000001", local_case_id="secret-case"))
        written = json.dumps(mod.dump(c, [mod.DUPLICATE_INITIALS, mod.INITIALS_CHECKED]))
        assert "260000001" not in written
        assert "secret-case" not in written


class TestDump:
    def test_absent_field_reports_zero(self):
        assert mod.dump(counters(), ["nope"])["nope"] == {"_total": 0, "_distinct": 0, "values": {}}

    def test_reports_every_value_most_common_first(self):
        c = counters()
        c["f"].update({"b": 3, "a": 5, "c": 1})
        out = mod.dump(c, ["f"])["f"]
        assert out["_total"] == 9
        assert out["_distinct"] == 3
        assert list(out["values"].items()) == [("a", 5), ("b", 3), ("c", 1)]


class TestAsUrl:
    def test_plain_postgresql_urls_are_routed_to_the_psycopg_dialect(self):
        # We ship psycopg 3 (grz-db's pin); a bare postgresql:// would make
        # SQLAlchemy look for psycopg2 instead.
        assert mod.as_url("postgresql://user@host/db") == "postgresql+psycopg://user@host/db"

    def test_an_explicit_driver_in_the_url_is_left_untouched(self):
        assert mod.as_url("postgresql+psycopg://u@h/db") == "postgresql+psycopg://u@h/db"
        assert mod.as_url("postgresql+psycopg2://u@h/db") == "postgresql+psycopg2://u@h/db"

    def test_non_postgresql_urls_pass_through(self):
        assert mod.as_url("sqlite:////tmp/x.sqlite") == "sqlite:////tmp/x.sqlite"

    def test_wraps_an_existing_path(self, tmp_path):
        p = tmp_path / "submission.db.sqlite"
        p.touch()
        assert mod.as_url(str(p)) == f"sqlite:///{p}"

    def test_exits_on_a_missing_path(self, tmp_path):
        with pytest.raises(SystemExit):
            mod.as_url(str(tmp_path / "nope.sqlite"))


class TestResolveDbUrl:
    def test_reads_a_grz_config_file(self, tmp_path):
        cfg = tmp_path / "config.yaml"
        cfg.write_text("db:\n  database_url: postgresql://user@host/grzdb\n")
        assert mod.resolve_db_url(config_file=str(cfg)) == "postgresql+psycopg://user@host/grzdb"

    def test_db_url_wins_over_a_config_file(self, tmp_path):
        assert mod.resolve_db_url("postgresql://user@host/db", "/nonexistent.yaml") == (
            "postgresql+psycopg://user@host/db"
        )

    def test_exits_when_the_config_has_no_url(self, tmp_path):
        cfg = tmp_path / "config.yaml"
        cfg.write_text("db: {}\n")
        with pytest.raises(SystemExit):
            mod.resolve_db_url(config_file=str(cfg))

    def test_exits_when_neither_is_given(self):
        with pytest.raises(SystemExit):
            mod.resolve_db_url()


class TestResolveGrzId:
    """The site id is read out of the submissions rather than retyped."""

    def test_taken_from_the_submissions(self):
        assert mod.resolve_grz_id(collections.Counter({"GRZK00123": 40})) == (
            "GRZK00123",
            mod.GRZ_ID_PATH,
        )

    def test_the_flag_overrides_what_the_data_says(self):
        found, source = mod.resolve_grz_id(collections.Counter({"GRZK00123": 40}), "GRZK00999")
        assert (found, source) == ("GRZK00999", "--grz-id")

    def test_exits_when_the_database_records_none(self):
        with pytest.raises(SystemExit) as exit:
            mod.resolve_grz_id(collections.Counter())
        assert "--grz-id" in str(exit.value)

    def test_exits_when_the_database_holds_more_than_one(self):
        with pytest.raises(SystemExit) as exit:
            mod.resolve_grz_id(collections.Counter({"GRZK00123": 40, "GRZK00456": 2}))
        message = str(exit.value)
        assert "more than one GRZ" in message
        # the counts are what tell a stray import from the site's own submissions
        assert "GRZK00123 (40)" in message and "GRZK00456 (2)" in message

    def test_the_flag_rescues_both_of_those(self):
        assert mod.resolve_grz_id(collections.Counter(), "GRZK00123")[0] == "GRZK00123"
        assert mod.resolve_grz_id(collections.Counter({"A": 1, "B": 1}), "GRZK00123")[0] == "GRZK00123"


class TestEndToEnd:
    def test_writes_a_report(self, sqlite_db, tmp_path):
        out = tmp_path / "report.json"
        mod.main(["--db-url", str(sqlite_db), "--out", str(out)])
        report = json.loads(out.read_text(encoding="utf-8"))

        assert report["grz_id"] == "GRZK00123"  # never passed on the command line
        assert report["grz_id_source"] == mod.GRZ_ID_PATH
        assert report["submissions_in_table"] == 4
        assert report["submissions_with_metadata"] == 2
        assert report["submissions_unparseable"] == 1  # the "{not json" row; the NULL row is not counted
        assert report["enum_fields"]["labData.libraryType"]["values"] == {"wgs": 1, "wxs": 1}
        assert report["freetext_fields"]["labData.labDataName"]["values"] == {"Blut DNA": 2}
        assert report["derived"]["tissueTypeId_is_BTO_format"]["values"] == {"yes": 2}
        # Two initial submissions, two different cases, one LE: nothing duplicated.
        assert report["derived"][mod.DUPLICATE_INITIALS]["values"] == {"0": 1}
        assert report["derived"][mod.INITIALS_CHECKED]["values"] == {
            "yes": 2,
            mod.NOT_QC_PASSED: 0,
            mod.REDACTED_CASE_ID: 0,
            mod.NOT_ATTRIBUTABLE: 0,
        }

    def test_the_duplicates_are_named_on_the_console_but_not_in_the_report(
        self, db_factory, submission, tmp_path, caplog
    ):
        # The console runs at the site that can go and look the case up; the
        # report is what leaves. Only one of the two may carry identifiers.
        db = db_factory([submission(submitter_id="260000001", local_case_id="ABC-17")] * 2)
        out = tmp_path / "report.json"
        with caplog.at_level("INFO"):
            mod.main(["--db-url", str(db), "--out", str(out)])
        assert "260000001: ABC-17 (2x)" in caplog.text
        assert "stay on this machine" in caplog.text
        written = out.read_text(encoding="utf-8")
        assert "260000001" not in written
        assert "ABC-17" not in written

    def test_a_clean_database_says_so(self, db_factory, submission, tmp_path, caplog):
        db = db_factory([submission(local_case_id="a"), submission(local_case_id="b")])
        with caplog.at_level("INFO"):
            mod.main(["--db-url", str(db), "--out", str(tmp_path / "report.json")])
        assert "no Leistungserbringer duplicated a QC-passed initial submission (1 checked)" in caplog.text

    def test_a_failed_submission_and_its_retry_are_not_a_duplicate(self, db_factory, submission, tmp_path):
        db = db_factory(
            [
                (submission(local_case_id="a"), False),  # rejected by QC
                (submission(local_case_id="a"), True),  # sent again, accepted
            ]
        )
        out = tmp_path / "report.json"
        mod.main(["--db-url", str(db), "--out", str(out)])
        report = json.loads(out.read_text(encoding="utf-8"))
        assert report["derived"][mod.DUPLICATE_INITIALS]["values"] == {"0": 1}
        assert report["derived"][mod.INITIALS_CHECKED]["values"][mod.NOT_QC_PASSED] == 1

    def test_the_case_id_comes_from_the_column_not_the_metadata(self, db_factory, submission, tmp_path):
        # Both documents carry the redaction placeholder; the pseudonym column
        # holds the real, distinct case ids. Reading the JSON would call these
        # two a duplicate pair.
        redacted = submission(local_case_id="REDACTED_LOCAL_CASE_ID")
        db = db_factory([(redacted, True, "real-1"), (redacted, True, "real-2")])
        out = tmp_path / "report.json"
        mod.main(["--db-url", str(db), "--out", str(out)])
        report = json.loads(out.read_text(encoding="utf-8"))
        assert report["derived"][mod.DUPLICATE_INITIALS]["values"] == {"0": 1}
        assert report["derived"][mod.INITIALS_CHECKED]["values"]["yes"] == 2

    def test_a_database_without_the_qc_column_says_so_and_counts_nothing(
        self, db_factory, submission, tmp_path, caplog
    ):
        # Reporting an unfiltered number as if it were filtered would be worse
        # than reporting none; every other count has to survive regardless.
        db = db_factory([submission(local_case_id="a")] * 2, with_qc_column=False)
        out = tmp_path / "report.json"
        with caplog.at_level("INFO"):
            mod.main(["--db-url", str(db), "--out", str(out)])
        assert mod.QC_COLUMN in caplog.text and mod.LOCAL_CASE_ID_COLUMN in caplog.text
        report = json.loads(out.read_text(encoding="utf-8"))
        assert report["derived"][mod.DUPLICATE_INITIALS]["values"] == {}
        assert report["derived"][mod.INITIALS_CHECKED]["values"][mod.NOT_QC_PASSED] == 2
        assert report["enum_fields"]["labData.libraryType"]["values"] == {"wgs": 2}

    def test_reports_duplicate_initial_submissions_per_le(self, db_factory, submission, tmp_path):
        db = db_factory(
            [
                submission(submitter_id="260000001", local_case_id="a"),
                submission(submitter_id="260000001", local_case_id="a"),  # the duplicate
                submission(submitter_id="260000001", local_case_id="b"),
                submission(submitter_id="260000002", local_case_id="a"),
            ]
        )
        out = tmp_path / "report.json"
        mod.main(["--db-url", str(db), "--out", str(out)])
        report = json.loads(out.read_text(encoding="utf-8"))
        # 260000001 duplicated one case, 260000002 duplicated nothing.
        assert report["derived"][mod.DUPLICATE_INITIALS]["values"] == {"1": 1, "0": 1}
        assert report["derived"][mod.INITIALS_CHECKED]["values"]["yes"] == 4
        assert "260000001" not in out.read_text(encoding="utf-8")

    def test_report_values_are_written_most_common_first(self, db_factory, submission, tmp_path):
        db = db_factory([submission(library_type="wxs"), submission(library_type="wxs"), submission()])
        out = tmp_path / "report.json"
        mod.main(["--db-url", str(db), "--out", str(out)])
        report = json.loads(out.read_text(encoding="utf-8"))
        # The order in the file is what the human review sees; it must be the
        # frequency order dump() builds, not alphabetical.
        assert list(report["enum_fields"]["labData.libraryType"]["values"]) == ["wxs", "wgs"]

    def test_default_filename_uses_the_derived_grz_id(self, sqlite_db, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        mod.main(["--db-url", str(sqlite_db)])
        assert list(tmp_path.glob("grz-survey-GRZK00123-*.json"))

    def test_the_flag_relabels_the_report_and_says_so(self, sqlite_db, tmp_path):
        out = tmp_path / "report.json"
        mod.main(["--db-url", str(sqlite_db), "--grz-id", "GRZK00999", "--out", str(out)])
        report = json.loads(out.read_text(encoding="utf-8"))
        assert report["grz_id"] == "GRZK00999"
        assert report["grz_id_source"] == "--grz-id"

    def test_a_database_without_the_id_needs_the_flag(self, db_factory, submission, tmp_path):
        db = db_factory([submission(grz_id=None)])
        with pytest.raises(SystemExit) as exit:
            mod.main(["--db-url", str(db)])
        assert mod.GRZ_ID_PATH in str(exit.value)

        out = tmp_path / "report.json"
        mod.main(["--db-url", str(db), "--grz-id", "GRZK00123", "--out", str(out)])
        assert json.loads(out.read_text(encoding="utf-8"))["grz_id"] == "GRZK00123"

    def test_a_database_mixing_two_grzs_is_refused(self, db_factory, submission):
        db = db_factory([submission(), submission(grz_id="GRZK00456")])
        with pytest.raises(SystemExit) as exit:
            mod.main(["--db-url", str(db)])
        assert "more than one GRZ" in str(exit.value)
