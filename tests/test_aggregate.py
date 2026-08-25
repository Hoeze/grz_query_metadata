from __future__ import annotations

import json

import pytest

from grz_query_metadata import aggregate as mod
from grz_query_metadata import ods


@pytest.mark.parametrize(
    "value,expected",
    [
        ("Illumina NovaSeq 6000", "illuminanovaseq6000"),
        ("illumina_novaseq-6000", "illuminanovaseq6000"),
        ("NovaSeq 6000 ", "novaseq6000"),
    ],
)
def test_normalise(value, expected):
    assert mod.normalise(value) == expected


def test_sheet_name_strips_illegal_characters_and_fits_the_31_character_cap():
    name = mod.sheet_name("donors[]/labData[]/sequenceData/bioinformaticsPipelineName")
    assert len(name) <= 31
    assert not set(name) & set("[]:*?/\\")
    assert name.endswith("bioinformaticsPipelineName")  # the tail is kept — it distinguishes labels


def test_row_order_sorts_numeric_labels_numerically():
    # The duplicate-per-LE distribution is keyed by a count, and its rows tie
    # constantly; alphabetically "10" would land between "1" and "2".
    assert sorted(["10", "2", "1"], key=mod.row_order) == ["1", "2", "10"]


def test_row_order_leaves_everything_else_alphabetical():
    assert sorted(["wgs", "1", "panel"], key=mod.row_order) == ["1", "panel", "wgs"]


def test_sheet_name_keeps_the_head_of_a_label_that_is_not_a_path():
    # A derived metric has no shared prefix to drop, so cutting its head leaves
    # gibberish: "duplicate_initial..." must not become "cate_initial...".
    assert mod.sheet_name("duplicate_initial_submissions_per_LE").startswith("duplicate_initial")
    assert mod.sheet_name("initial_submissions_checked_for_duplicates").startswith("initial_sub")


def test_sheet_name_leaves_a_short_label_alone():
    assert mod.sheet_name("labData.libraryType") == "labData.libraryType"


def test_sheet_name_never_comes_back_empty():
    assert mod.sheet_name("[]:*?/\\") == "field"


def report(label, **fields):
    """A survey report as it looks in memory, i.e. after load_reports has
    tagged it with the column heading it gets."""
    return {
        "script_version": "1.1.0",
        "grz_id": label,
        "grz_id_source": "submission/genomicDataCenterId",
        "_label": label,
        "generated": "2026-08-13",
        "submissions_in_table": 10,
        "submissions_with_metadata": 9,
        "submissions_unparseable": 0,
        "duplicate_check_available": True,
        "enum_fields": {},
        "freetext_fields": {},
        "derived": {},
        **fields,
    }


def write(tmp_path, name, payload):
    """The same report on disk, where the label does not exist yet."""
    p = tmp_path / name
    p.write_text(json.dumps({k: v for k, v in payload.items() if k != "_label"}), encoding="utf-8")
    return str(p)


def field(values):
    return {"_total": sum(values.values()), "_distinct": len(values), "values": values}


def values_of(rows):
    return [list(r.values) for r in rows]


class TestField:
    """The per-field model both the sheet and its index row are derived from."""

    def make(self, per_grz, **kwargs):
        return mod.Field(section="enum_fields", label="labData.libraryType", per_grz=per_grz, **kwargs)

    def test_totals_add_up_across_grzs(self):
        f = self.make({"A": {"wgs": 3, "wes": 1}, "B": {"wgs": 2}})
        assert f.total == 6
        assert f.total_for("wgs") == 5
        assert f.observed == {"wgs", "wes"}

    def test_a_grz_that_never_reported_the_field_is_still_a_column(self):
        f = self.make({"A": {"wgs": 3}, "B": {}})
        assert f.grzs == ["A", "B"]
        assert f.total_for("wgs") == 3

    def test_values_are_ordered_by_frequency_then_alphabetically(self):
        f = self.make({"A": {"wes": 1, "wgs": 1, "panel": 5}})
        assert f.values == ["panel", "wes", "wgs"]

    def test_declared_values_nobody_used_still_get_a_row_at_the_end(self):
        f = self.make({"A": {"wgs": 3}}, declared=["wgs", "wes", "panel"])
        assert f.values == ["wgs", "panel", "wes"]  # the unused two sort alphabetically
        assert f.unused == ["wes", "panel"]  # ... but are listed in schema order

    def test_unused_is_empty_without_a_schema(self):
        assert self.make({"A": {"wgs": 3}}).unused == []

    def test_vocabulary_coverage_counts_observations_not_values(self):
        f = self.make({"A": {"wgs": 9, "wes": 1}}, vocab={"wgs"})
        assert f.covered == 9  # one value of two, but nine observations of ten
        assert f.in_vocab("wgs") and not f.in_vocab("wes")


class TestSheetFunctions:
    """Each sheet is built by one function, so each can be read on its own."""

    def test_summary_row_per_report_then_a_total(self):
        rows = mod.summary_rows([report("GRZ_A"), report("GRZ_B")])
        source = "submission/genomicDataCenterId"
        assert values_of(rows)[1:] == [
            ["GRZ_A", 10, 9, 0, "yes", source, "1.1.0", "2026-08-13"],
            ["GRZ_B", 10, 9, 0, "yes", source, "1.1.0", "2026-08-13"],
            [],
            ["TOTAL", 20, 18, 0],
        ]

    def test_summary_warns_only_when_versions_differ(self):
        same = mod.summary_rows([report("GRZ_A"), report("GRZ_B")])
        assert not any("WARNING" in str(r.values) for r in same)

        mixed = mod.summary_rows([report("GRZ_A"), report("GRZ_B", script_version="2.0.0")])
        assert mixed[-1].style == ods.PROBLEM
        assert "1.1.0" in mixed[-1].values[0] and "2.0.0" in mixed[-1].values[0]

    @pytest.mark.parametrize(
        "flag,expected",
        [
            ({"duplicate_check_available": True}, "yes"),
            ({"duplicate_check_available": False}, "no"),
            ({}, "unknown"),
        ],
    )
    def test_summary_says_whether_the_duplicate_check_could_run(self, flag, expected):
        # A report from before the check existed carries no such key at all.
        r = report("GRZ_A", **flag)
        if not flag:
            del r["duplicate_check_available"]
        assert mod.duplicate_check_status(r) == expected
        assert values_of(mod.summary_rows([r]))[1][4] == expected

    def test_summary_warns_when_a_grz_never_ran_the_duplicate_check(self):
        # Its sheet shows zeroes either way; only this row tells them apart.
        clean = mod.summary_rows([report("GRZ_A"), report("GRZ_B")])
        assert not any("duplicate initial submissions" in str(r.values) for r in clean)

        mixed = mod.summary_rows([report("GRZ_A"), report("GRZ_B", duplicate_check_available=False)])
        assert mixed[-1].style == ods.PROBLEM
        warning = mixed[-1].values[0]
        assert "GRZ_B" in warning and "GRZ_A" not in warning
        assert "not checked" in warning

    def test_index_leaves_the_schema_columns_empty_without_one(self):
        fields = mod.collect_fields([report("A", enum_fields={"labData.libraryType": field({"wgs": 1})})])
        assert values_of(mod.index_rows(fields))[1] == [
            "enum_fields",
            "labData.libraryType",
            1,
            "",
            "",
            "labData.libraryType",
        ]

    def test_index_highlights_fields_with_unused_declared_values(self):
        fields = mod.collect_fields(
            [report("A", enum_fields={"labData.libraryType": field({"wgs": 1})})],
            declared={"labData.libraryType": ["wgs", "wes"]},
        )
        (row,) = mod.index_rows(fields)[1:]
        assert list(row.values) == ["enum_fields", "labData.libraryType", 1, 2, 1, "labData.libraryType"]
        assert row.styles == {4: ods.UNUSED}

    def test_field_header_widens_only_with_the_matching_option(self):
        per_grz = {"A": {"wgs": 1}}
        plain = mod.Field("enum_fields", "f", per_grz)
        assert mod.field_header(plain) == ["value", "A", "TOTAL", "share"]
        assert mod.field_header(mod.Field("enum_fields", "f", per_grz, declared=["wgs"]))[-1] == (
            "declared in schema"
        )
        assert mod.field_header(mod.Field("enum_fields", "f", per_grz, vocab={"wgs"}))[-1] == (
            "in proposed vocabulary"
        )

    def test_field_rows_end_with_the_footers_the_options_add(self):
        f = mod.Field("enum_fields", "f", {"A": {"wgs": 3}}, declared=["wgs", "wes"], vocab={"wgs"})
        rows = values_of(mod.field_rows(f))
        assert rows[0] == ["value", "A", "TOTAL", "share", "in proposed vocabulary", "declared in schema"]
        assert rows[1] == ["wgs", 3, 3, 1.0, "yes", "yes"]
        assert rows[2] == ["wes", 0, 0, 0, "NO", "yes"]
        assert ["TOTAL", 3, 3] in rows
        assert ["declared but never used (1 of 2)", "wes"] in rows
        assert ["coverage by proposed vocabulary", "3/3 observations", "100.0%"] in rows

    def test_a_declared_but_unused_row_is_flagged_as_unused_not_as_off_vocabulary(self):
        f = mod.Field("enum_fields", "f", {"A": {"wgs": 3}}, declared=["wgs", "wes"], vocab={"wgs"})
        assert mod.field_rows(f)[2].style == ods.UNUSED


class TestSpreadsheet:
    def test_one_column_per_report_and_a_total(self, tmp_path, read_ods):
        a = write(
            tmp_path,
            "a.json",
            report("GRZ_A", enum_fields={"labData.libraryType": field({"wgs": 3, "wes": 1})}),
        )
        b = write(
            tmp_path,
            "b.json",
            report("GRZ_B", enum_fields={"labData.libraryType": field({"wgs": 2})}),
        )
        out = tmp_path / "usage.ods"
        mod.main([a, b, "-o", str(out)])

        rows = read_ods(out)["labData.libraryType"]
        assert rows[0] == ["value", "GRZ_A", "GRZ_B", "TOTAL", "share"]
        assert rows[1] == ["wgs", 3, 2, 5, round(5 / 6, 4)]
        assert rows[2] == ["wes", 1, 0, 1, round(1 / 6, 4)]
        assert rows[-1] == ["TOTAL", 4, 2, 6]  # per-GRZ totals, then all observations

    def test_summary_and_index_sheets_are_present(self, tmp_path, read_ods):
        p = write(
            tmp_path,
            "a.json",
            report("GRZ_A", enum_fields={"labData.libraryType": field({"wgs": 3})}),
        )
        out = tmp_path / "usage.ods"
        mod.main([p, "-o", str(out)])

        sheets = read_ods(out)
        assert list(sheets) == ["summary", "index", "labData.libraryType"]
        assert sheets["summary"][1] == [
            "GRZ_A",
            10,
            9,
            0,
            "yes",
            "submission/genomicDataCenterId",
            "1.1.0",
            "2026-08-13",
        ]
        # the declared/never-used columns stay empty without --schema
        assert sheets["index"][1] == [
            "enum_fields",
            "labData.libraryType",
            1,
            None,
            None,
            "labData.libraryType",
        ]

    def test_label_falls_back_to_the_filename_without_a_grz_id(self, tmp_path, read_ods):
        p = write(tmp_path, "grz-survey-2026-08-13.json", report(None))
        out = tmp_path / "usage.ods"
        mod.main([p, "-o", str(out)])
        assert read_ods(out)["summary"][1][0] == "grz-survey-2026-08-13"

    def test_mixed_script_versions_are_flagged(self, tmp_path, read_ods, cell_styles):
        a = write(tmp_path, "a.json", report("GRZ_A"))
        b = write(tmp_path, "b.json", report("GRZ_B", script_version="2.0.0"))
        out = tmp_path / "usage.ods"
        mod.main([a, b, "-o", str(out)])

        assert "WARNING" in read_ods(out)["summary"][-1][0]
        assert cell_styles(out, "summary")[-1] == [ods.PROBLEM]

    def test_schema_adds_a_zero_row_for_values_nobody_used(self, tmp_path, schema_dir, read_ods, cell_styles):
        p = write(
            tmp_path,
            "a.json",
            report("GRZ_A", enum_fields={"labData.libraryType": field({"wgs": 4})}),
        )
        out = tmp_path / "usage.ods"
        mod.main([p, "-o", str(out), "--schema", str(schema_dir / "grz-schema.json")])

        rows = read_ods(out)["labData.libraryType"]
        assert rows[1] == ["wgs", 4, 4, 1.0, "yes"]
        # "wes" is declared but was never submitted: a zero row, still marked declared
        assert rows[2] == ["wes", 0, 0, 0, "yes"]
        assert cell_styles(out, "labData.libraryType")[2] == [ods.UNUSED] * 5
        assert rows[-1] == ["declared but never used (1 of 2)", "wes"]

    def test_vocab_flags_values_outside_the_proposed_enum(self, tmp_path, read_ods, cell_styles):
        vocab = tmp_path / "instrument-model.json"
        vocab.write_text(json.dumps({"enum": ["Illumina NovaSeq 6000"]}), encoding="utf-8")
        p = write(
            tmp_path,
            "a.json",
            report(
                "GRZ_A",
                freetext_fields={
                    "labData.sequencerModel": field({"illumina novaseq-6000": 3, "PacBio Revio": 1})
                },
            ),
        )
        out = tmp_path / "usage.ods"
        mod.main([p, "-o", str(out), f"--vocab=labData.sequencerModel={vocab}"])

        rows = read_ods(out)["labData.sequencerModel"]
        assert rows[0][-1] == "in proposed vocabulary"
        assert rows[1][-1] == "yes"  # matched loosely, despite the different spelling
        assert rows[2][-1] == "NO"
        assert cell_styles(out, "labData.sequencerModel")[2] == [ods.OUTSIDE_VOCAB] * 5
        assert rows[-1] == ["coverage by proposed vocabulary", "3/4 observations", "75.0%"]

    def test_an_empty_vocabulary_still_produces_the_coverage_column(self, tmp_path, read_ods):
        # An empty enum is a real input — it covers nothing — and must show
        # 0% rather than silently switching the analysis off.
        vocab = tmp_path / "empty.json"
        vocab.write_text(json.dumps({"enum": []}), encoding="utf-8")
        p = write(
            tmp_path,
            "a.json",
            report("GRZ_A", freetext_fields={"labData.sequencerModel": field({"NovaSeq": 3})}),
        )
        out = tmp_path / "usage.ods"
        mod.main([p, "-o", str(out), f"--vocab=labData.sequencerModel={vocab}"])

        rows = read_ods(out)["labData.sequencerModel"]
        assert rows[0][-1] == "in proposed vocabulary"
        assert rows[1][-1] == "NO"
        assert rows[-1] == ["coverage by proposed vocabulary", "0/3 observations", "0.0%"]

    def test_a_vocab_for_an_unknown_field_is_refused(self, tmp_path):
        # A typo'd field label must not silently skip the check while the
        # console claims coverage was measured.
        vocab = tmp_path / "v.json"
        vocab.write_text(json.dumps({"enum": ["x"]}), encoding="utf-8")
        p = write(tmp_path, "a.json", report("GRZ_A"))
        with pytest.raises(SystemExit) as exit:
            mod.main([p, "-o", str(tmp_path / "usage.ods"), f"--vocab=labdata.sequencerModel={vocab}"])
        assert "labdata.sequencerModel" in str(exit.value)

    def test_a_vocab_spec_without_the_field_prefix_is_refused(self):
        with pytest.raises(SystemExit) as exit:
            mod.load_vocabs(["instrument-model.json"])
        assert "FIELD=" in str(exit.value)

    def test_a_vocab_file_without_an_enum_is_refused(self, tmp_path):
        bad = tmp_path / "bad.json"
        bad.write_text(json.dumps({"values": []}), encoding="utf-8")
        with pytest.raises(SystemExit) as exit:
            mod.load_vocabs([f"f={bad}"])
        assert '"enum"' in str(exit.value)

    def test_two_reports_with_the_same_label_are_refused(self, tmp_path):
        # Same grz_id twice would merge into one column while the summary
        # counted both — the spreadsheet would contradict itself.
        a = write(tmp_path, "a.json", report("GRZ_A"))
        b = write(tmp_path, "b.json", report("GRZ_A"))
        with pytest.raises(SystemExit) as exit:
            mod.main([a, b, "-o", str(tmp_path / "usage.ods")])
        assert "GRZ_A" in str(exit.value)

    def test_colliding_sheet_names_are_disambiguated(self, tmp_path, read_ods):
        # Both labels truncate to the same last-31 characters; two same-named
        # tables would make the .ods invalid and hide one field.
        long_a, long_b = "a" + "x" * 34, "b" + "x" * 34
        p = write(
            tmp_path,
            "a.json",
            report("GRZ_A", enum_fields={long_a: field({"v": 1}), long_b: field({"w": 1})}),
        )
        out = tmp_path / "usage.ods"
        mod.main([p, "-o", str(out)])

        sheets = read_ods(out)
        names = [n for n in sheets if n not in ("summary", "index")]
        assert len(names) == 2 and len(set(names)) == 2
        assert {r[5] for r in sheets["index"][1:]} == set(names)  # index points at real sheets


class TestSchemaLoading:
    def test_a_dangling_ref_degrades_to_a_warning_not_a_crash(self, tmp_path, caplog):
        schema = tmp_path / "grz-schema.json"
        schema.write_text(
            json.dumps({"properties": {"submission": {"$ref": "#/$defs/Gone"}}}), encoding="utf-8"
        )
        assert mod.load_declared(schema) == {}
        assert "could not walk the schema" in caplog.text

    def test_schema_enums_the_survey_misses_are_warned_about(self, tmp_path, caplog):
        schema = tmp_path / "grz-schema.json"
        schema.write_text(
            json.dumps(
                {
                    "properties": {
                        "submission": {
                            "properties": {
                                "submissionType": {"enum": ["initial"]},  # surveyed
                                "brandNewField": {"enum": ["a", "b"]},  # not surveyed
                            }
                        }
                    }
                }
            ),
            encoding="utf-8",
        )
        mod.load_declared(schema)
        assert "does not cover" in caplog.text
        assert "submission/brandNewField" in caplog.text
