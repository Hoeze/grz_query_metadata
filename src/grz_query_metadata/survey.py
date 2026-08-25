"""Survey controlled-vocabulary usage in GRZ submission metadata.

Reads the `submissions` table of a GRZ internal submission database and counts
how often each value of each surveyed field occurs. Works against SQLite and
PostgreSQL alike: the metadata JSON is selected as a whole and inspected in
Python, so no dialect-specific JSON SQL is involved.

Emits a single JSON report. The report contains counts and, for technical
fields, the distinct values found. It contains no tanG, no donor pseudonym, no
localCaseId, no submitterId, no file paths and no dates. submitterId and
localCaseId are read, but only to group rows in memory for the duplicate
initial submission check; see :class:`InitialSubmissions`. The localCaseId
comes from the `pseudonym` column, which is where grz-db keeps it, not from a
donor pseudonym; see :data:`LOCAL_CASE_ID_COLUMN`.

Usage
-----
    grz-survey-metadata --db-url sqlite:////path/to/submission.db.sqlite
    grz-survey-metadata --db-url postgresql://user@host/grzdb

    # read the URL from a grz config file instead
    grz-survey-metadata --config-file /etc/grz/config.yaml
"""

from __future__ import annotations

import argparse
import collections
import datetime
import json
import logging
import os
import sys
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

import yaml
from sqlalchemy import create_engine, inspect, text
from sqlalchemy.exc import SQLAlchemyError

from . import __version__, setup_cli_logging
from .fields import BTO_ID, ENUM_FIELDS, FREETEXT_FIELDS, split_segment

log = logging.getLogger(__name__)

Counters = dict[str, collections.Counter[str]]

# Every submission records the GRZ that received it, so the report can label
# itself instead of asking for the id to be retyped. This is the site's own
# identifier, not anyone's patient or institution.
GRZ_ID_PATH = "submission/genomicDataCenterId"

# Read to group rows in memory, never counted and never written to the report:
# the duplicate check has to know which LE sent a submission and which case it
# was for, and both of those identify a real institution or a real episode of
# care. submitterId is the §293 SGB V id of the submitting institution, and it
# is not one of the fields grz-db redacts, so the metadata document still has it.
SUBMITTER_ID_PATH = "submission/submitterId"

# Taken from the surveyed field rather than retyped, so the duplicate check and
# the submissionType counts cannot end up reading different places.
SUBMISSION_TYPE_PATH = ENUM_FIELDS["submission.submissionType"]

# grz-db's own definition of a QC-passed initial submission, which its
# _qc_passed_initial_of() and its one-initial-per-case index both use. Only
# basic QC: detailed QC runs on a selected sample, so filtering on it would
# throw away most submissions rather than the ones that failed.
QC_COLUMN = "basic_qc_passed"

# The submitter's localCaseId, which grz-db stores in a column named for the
# other redacted field rather than for itself (its _METADATA_FIELD_TO_COLUMN
# maps local_case_id -> pseudonym). It is NOT the donor pseudonym: that lives in
# the separate donors table and is never touched here.
#
# It has to come from the column, because grz-db redacts localCaseId inside the
# stored metadata document. Reading the JSON would key every redacted submission
# of an LE on one placeholder and report them all as duplicates of each other.
LOCAL_CASE_ID_COLUMN = "pseudonym"

# grz-pydantic-models' LOCAL_CASE_ID_PLACEHOLDERS. A placeholder identifies no
# case, so it must never be keyed on, which is also how grz-db's
# SubmitterLocalCaseResolver treats it.
LOCAL_CASE_ID_PLACEHOLDERS = ("", "REDACTED_LOCAL_CASE_ID")

DUPLICATE_INITIALS = "duplicate_initial_submissions_per_LE"
INITIALS_CHECKED = "initial_submissions_checked_for_duplicates"

# Why an initial submission was left out. Written even at zero, so that the
# count above can be read as a total rather than a floor.
NOT_QC_PASSED = "no - basic QC not passed"
REDACTED_CASE_ID = "no - localCaseId redacted"
NOT_ATTRIBUTABLE = "no - no submitterId or localCaseId"


# --------------------------------------------------------------------------
# JSON path walking
# --------------------------------------------------------------------------


def walk(node: Any, parts: list[str]):
    """Yield every value reachable by `parts`, descending into [] segments."""
    if node is None:
        return
    if not parts:
        yield node
        return
    head, rest = parts[0], parts[1:]
    key, is_array = split_segment(head)
    if is_array:
        seq = node.get(key) if isinstance(node, dict) else None
        if isinstance(seq, list):
            for item in seq:
                yield from walk(item, rest)
    else:
        if isinstance(node, dict) and key in node:
            yield from walk(node[key], rest)


def collect(meta: dict, path: str) -> list[Any]:
    return [v for v in walk(meta, path.split("/")) if v is not None]


def first(meta: dict, path: str) -> str | None:
    """The value at a single-valued path, as a string, or None if absent.

    Goes through collect() so that a malformed document yields None instead of
    raising, exactly as it does for the counted fields.
    """
    found = collect(meta, path)
    return str(found[0]) if found else None


def as_bool(value: Any) -> bool | None:
    """A tri-state column value as a bool, NULL staying None.

    SQLite has no boolean type and a text() query carries no type information
    for SQLAlchemy to convert with, so the driver hands back 1/0 there and
    True/False on PostgreSQL. `is True` has to mean the same thing on both.
    """
    return None if value is None else bool(value)


# --------------------------------------------------------------------------
# Derived checks: things a plain value count cannot answer
# --------------------------------------------------------------------------


def derived_metrics(meta: dict, d: Counters) -> None:
    # --- tissueTypeId: does it look like a BTO identifier? ---
    # collect() rather than a hand-rolled walk: it shrugs off malformed shapes
    # (a null donor, labData not a list) exactly like the field counting does,
    # and it reads the path from the same table, so the two cannot diverge.
    for tid in collect(meta, FREETEXT_FIELDS["labData.tissueTypeId"]):
        d["tissueTypeId_is_BTO_format"]["yes" if BTO_ID.match(str(tid)) else "no"] += 1


class InitialSubmissions:
    """Which Leistungserbringer sent a QC-passed initial submission for which case.

    A case is meant to arrive once as `initial`; what follows it is a follow-up
    or a correction. A second `initial` for the same
    (submitterId, localCaseId) is therefore a re-submission nobody intended,
    and how many of those each LE produced is the metric.

    Only submissions that passed basic QC count. One that failed was rejected
    and the LE was meant to send it again, so counting the replacement as a
    duplicate would report the process working as a fault. grz-db draws the
    line in the same place: its one-initial-per-case index constrains QC-passed
    initials only, which is why this check finds anything at all — what it
    surfaces is the part that index does not cover, namely submissions from
    before it existed and those with no case linked.

    This is the one check that cannot be a per-submission
    :func:`derived_metrics` entry: duplication only becomes visible once the
    whole table has been read.

    The case id comes from the `pseudonym` column rather than the metadata
    document, because grz-db redacts localCaseId in the JSON it stores; see
    :data:`LOCAL_CASE_ID_COLUMN`.

    Neither identifier leaves this object. What reaches the report is the
    distribution alone — how many LEs duplicated nothing, how many duplicated
    once, and so on — keyed by the exact count and never by the LE. Every
    initial submission left out is counted under the reason it was left out,
    rather than quietly lowering the duplicate figure.
    """

    def __init__(self) -> None:
        self.cases: collections.Counter[tuple[str, str]] = collections.Counter()
        self.skipped: collections.Counter[str] = collections.Counter()

    def add(self, meta: dict, basic_qc_passed: bool | None, local_case_id: str | None) -> None:
        """Fold one submission in. Anything but an initial submission is ignored.

        `basic_qc_passed` and `local_case_id` are the columns, not the metadata
        document: None means the row holds no value, or that this database has
        no such column. For QC that is not a pass, so it is not counted. For the
        case id it is nothing to group by, so it is not counted either.
        """
        if first(meta, SUBMISSION_TYPE_PATH) != "initial":
            return
        if basic_qc_passed is not True:
            self.skipped[NOT_QC_PASSED] += 1
            return
        if local_case_id in LOCAL_CASE_ID_PLACEHOLDERS:
            self.skipped[REDACTED_CASE_ID] += 1
            return
        submitter = first(meta, SUBMITTER_ID_PATH)
        if submitter is None or local_case_id is None:
            self.skipped[NOT_ATTRIBUTABLE] += 1
            return
        self.cases[(submitter, local_case_id)] += 1

    def duplicates_per_le(self) -> collections.Counter[str]:
        """LE -> initial submissions it sent beyond the first, summed over its cases.

        Every LE that sent anything at all appears, zeroes included: "most LEs
        duplicate nothing" is half of what the distribution has to say.
        """
        per_le: collections.Counter[str] = collections.Counter()
        for (submitter, _case), seen in self.cases.items():
            per_le[submitter] += seen - 1
        return per_le

    def duplicated_cases(self) -> dict[str, list[tuple[str, int]]]:
        """LE -> its duplicated cases and how many times each was sent, worst first.

        The identifiers themselves, which is what :meth:`write` deliberately
        strips. Only the console ever sees this: it runs at the site that can
        go and look the cases up, and it never reaches the report.
        """
        by_le: dict[str, list[tuple[str, int]]] = collections.defaultdict(list)
        for (submitter, case), seen in self.cases.items():
            if seen > 1:
                by_le[submitter].append((case, seen))
        for cases in by_le.values():
            cases.sort(key=lambda pair: (-pair[1], pair[0]))
        return dict(sorted(by_le.items(), key=lambda item: (-sum(n for _, n in item[1]), item[0])))

    def write(self, counters: Counters) -> None:
        """Add the two derived counters to `counters`, identifiers left behind."""
        counters[DUPLICATE_INITIALS].update(str(n) for n in self.duplicates_per_le().values())
        checked = counters[INITIALS_CHECKED]
        checked["yes"] += sum(self.cases.values())
        # Explicit zeroes: a reason sitting at 0 is what says the count above
        # was computed from every initial submission rather than from whichever
        # part happened to be usable.
        for reason in (NOT_QC_PASSED, REDACTED_CASE_ID, NOT_ATTRIBUTABLE):
            checked[reason] += self.skipped[reason]


# --------------------------------------------------------------------------
# Counting
# --------------------------------------------------------------------------


def count_submission(meta: dict, counters: Counters) -> None:
    """Fold one submission's metadata document into `counters`."""
    for label, path in (*ENUM_FIELDS.items(), *FREETEXT_FIELDS.items()):
        for v in collect(meta, path):
            counters[label][str(v)] += 1
    derived_metrics(meta, counters)


def dump(counters: Counters, labels: Iterable[str]) -> dict[str, dict]:
    """Render the counters for `labels`, most common value first."""
    out: dict[str, dict] = {}
    for label in labels:
        c = counters.get(label)
        if not c:
            out[label] = {"_total": 0, "_distinct": 0, "values": {}}
            continue
        out[label] = {
            "_total": sum(c.values()),
            "_distinct": len(c),
            "values": dict(c.most_common()),
        }
    return out


@dataclass
class SurveyResult:
    """What one pass over a submissions table yields."""

    counters: Counters
    n_rows: int  # rows in the table
    n_with_metadata: int  # rows carrying a usable metadata document
    n_unparseable: int  # rows whose metadata could not be parsed
    grz_ids: collections.Counter[str]  # genomicDataCenterId values seen, with counts
    initials: InitialSubmissions  # duplicate check state; the report gets only its counts
    duplicate_check_available: bool  # whether the table has the columns the duplicate check needs


DUPLICATE_CHECK_COLUMNS = (QC_COLUMN, LOCAL_CASE_ID_COLUMN)


def missing_columns(engine, wanted: Sequence[str]) -> list[str]:
    """Which of `wanted` the submissions table does not have.

    grz-db grew these columns over time, so a database old enough to predate one
    must degrade to "the duplicate check cannot run" — checked up front rather
    than by letting the SELECT fail, since on PostgreSQL a failed statement
    aborts the transaction the rest of the pass needs.
    """
    try:
        present = {c["name"] for c in inspect(engine).get_columns("submissions")}
    except SQLAlchemyError as e:
        log.warning("could not inspect the columns of the submissions table: %r", e)
        return list(wanted)
    return [c for c in wanted if c not in present]


def survey(db_url: str) -> SurveyResult:
    """Count every surveyed field across the `submissions` table."""
    engine = create_engine(db_url)
    counters: Counters = collections.defaultdict(collections.Counter)
    grz_ids: collections.Counter[str] = collections.Counter()
    initials = InitialSubmissions()
    n_rows = n_with_metadata = n_unparseable = 0

    absent = missing_columns(engine, DUPLICATE_CHECK_COLUMNS)
    can_check = not absent
    if absent:
        log.warning(
            "this database has no %s column, so the duplicate initial submission check "
            "reports nothing. Every other count is unaffected.",
            " and no ".join(absent),
        )

    with engine.connect() as conn:
        # No JSON SQL: the driver hands back dict (JSON/JSONB) or str (older rows).
        # The two extra columns are the only thing read outside the JSON, and
        # they are selected only when the duplicate check can use both.
        extra = "".join(f", {c}" for c in DUPLICATE_CHECK_COLUMNS) if can_check else ""
        result = conn.execute(text(f"SELECT id, submission_metadata{extra} FROM submissions"))
        for row in result.mappings():
            n_rows += 1
            meta = row["submission_metadata"]
            if meta is None:
                continue
            if isinstance(meta, str):
                try:
                    meta = json.loads(meta)
                except json.JSONDecodeError:
                    n_unparseable += 1
                    continue
            if not isinstance(meta, dict):
                n_unparseable += 1
                continue
            n_with_metadata += 1
            count_submission(meta, counters)
            # Separate from count_submission because it is not a fold into a
            # counter: it keeps state of its own until the table is exhausted,
            # and it is the one check that reads a column rather than the JSON.
            initials.add(
                meta,
                as_bool(row[QC_COLUMN]) if can_check else None,
                row[LOCAL_CASE_ID_COLUMN] if can_check else None,
            )
            grz_ids.update(str(v) for v in collect(meta, GRZ_ID_PATH))

    initials.write(counters)
    return SurveyResult(counters, n_rows, n_with_metadata, n_unparseable, grz_ids, initials, can_check)


def resolve_grz_id(seen: collections.Counter[str], override: str | None = None) -> tuple[str, str]:
    """The site id for the report, and where it came from.

    The submissions record it themselves, so it is read out of the data rather
    than retyped. `--grz-id` overrides that, and is needed when the database
    holds no id at all or — which would be worth investigating — the
    submissions of more than one GRZ.
    """
    if override:
        return override, "--grz-id"
    if not seen:
        sys.exit(
            f"no {GRZ_ID_PATH} found in any submission of this database.\n"
            f"Pass --grz-id GRZKxxxxx to label the report yourself."
        )
    if len(seen) > 1:
        found = ", ".join(f"{gid} ({n})" for gid, n in seen.most_common())
        sys.exit(
            f"this database holds submissions from more than one GRZ: {found}.\n"
            f"Pass --grz-id GRZKxxxxx to say which one this report is for."
        )
    return next(iter(seen)), GRZ_ID_PATH


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def as_url(value: str) -> str:
    """Accept either a SQLAlchemy URL or a plain path to a SQLite file."""
    if "://" in value:
        # Plain postgresql:// would make SQLAlchemy look for psycopg2; we ship
        # psycopg 3 (the driver grz-db pins), so route to its dialect. An
        # explicit +driver in the URL is left untouched.
        if value.startswith("postgresql://"):
            return "postgresql+psycopg://" + value.removeprefix("postgresql://")
        return value
    path = os.path.abspath(os.path.expanduser(value))
    if not os.path.exists(path):
        sys.exit(
            f"no such file: {path}\n"
            f"(pass a SQLAlchemy URL such as postgresql://user@host/db "
            f"if this is not a SQLite file)"
        )
    return f"sqlite:///{path}"


def resolve_db_url(db_url: str | None = None, config_file: str | None = None) -> str:
    """The database URL, taken from `db_url` or from a grz config file."""
    if db_url:
        return as_url(db_url)
    if not config_file:
        sys.exit("either db_url or config_file is required")
    with open(config_file, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    url = (cfg.get("db") or {}).get("database_url")
    if not url:
        sys.exit(f"no db.database_url found in {config_file}")
    return as_url(url)


def build_report(result: SurveyResult, grz_id: str, grz_id_source: str = GRZ_ID_PATH) -> dict:
    """The JSON document a GRZ sends back."""
    counters = result.counters
    return {
        "script_version": __version__,
        "grz_id": grz_id,
        # Read out of the submissions, or overridden on the command line. Worth
        # recording, since only one of the two can have been mistyped.
        "grz_id_source": grz_id_source,
        "generated": datetime.date.today().isoformat(),
        "submissions_in_table": result.n_rows,
        "submissions_with_metadata": result.n_with_metadata,
        "submissions_unparseable": result.n_unparseable,
        # Whether this database has the columns the duplicate check reads.
        # Without them it counts nothing, and the aggregation has no other way
        # to tell that apart from a site that simply has no duplicates.
        "duplicate_check_available": result.duplicate_check_available,
        "enum_fields": dump(counters, ENUM_FIELDS),
        "freetext_fields": dump(counters, FREETEXT_FIELDS),
        "derived": dump(counters, sorted(set(counters) - set(ENUM_FIELDS) - set(FREETEXT_FIELDS))),
    }


def report_filename(generated: str, grz_id: str) -> str:
    return f"grz-survey-{grz_id}-{generated}.json"


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="grz-survey-metadata",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--db-url", help="SQLAlchemy URL, e.g. sqlite:///... or postgresql://...")
    src.add_argument("--config-file", help="grz config file containing db.database_url")
    ap.add_argument(
        "--grz-id",
        default=None,
        help=f"your site id, e.g. GRZK00123. Only needed to override what the "
        f"submissions themselves record in {GRZ_ID_PATH}, or to supply it when "
        f"they record nothing.",
    )
    ap.add_argument("--out", default=None, help="output file (default: grz-survey-<grzid>-<date>.json)")
    return ap


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    setup_cli_logging()

    db_url = resolve_db_url(args.db_url, args.config_file)
    result = survey(db_url)
    grz_id, source = resolve_grz_id(result.grz_ids, args.grz_id)
    report = build_report(result, grz_id, source)

    out = args.out or report_filename(report["generated"], grz_id)
    with open(out, "w", encoding="utf-8") as fh:
        # No sort_keys: it would alphabetise the values sections and destroy
        # the most-common-first ordering that makes the human review possible.
        # The report is built deterministically, so the file diffs fine as is.
        json.dump(report, fh, indent=2, ensure_ascii=False)

    log.info("reporting as %s (from %s)", grz_id, source)
    log.info("%d of %d submissions had metadata; wrote %s", result.n_with_metadata, result.n_rows, out)
    if result.n_unparseable:
        log.warning("%d rows had metadata that could not be parsed", result.n_unparseable)

    # The duplicate finding, said out loud. The spreadsheet is where it gets
    # compared across GRZs, but the site that can act on it is this one — and
    # this is the only place the identifiers may appear, so it names them.
    duplicated = result.initials.duplicated_cases()
    n_les = len(result.initials.duplicates_per_le())
    if duplicated:
        detail = "\n".join(
            f"    {submitter}: " + ", ".join(f"{case} ({n}x)" for case, n in cases)
            for submitter, cases in duplicated.items()
        )
        log.warning(
            "%d of %d Leistungserbringer sent the same case as a QC-passed initial "
            "submission more than once:\n"
            "  (these submitterIds and localCaseIds stay on this machine — the report "
            "carries only how many duplicates each LE had, and no identifier)\n%s",
            len(duplicated),
            n_les,
            detail,
        )
    elif n_les:
        log.info("no Leistungserbringer duplicated a QC-passed initial submission (%d checked)", n_les)
    unattributable = result.initials.skipped[NOT_ATTRIBUTABLE]
    if unattributable:
        log.warning(
            "%d QC-passed initial submission(s) carried no submitterId or no localCaseId, "
            "so they could not be checked for duplicates",
            unattributable,
        )
    log.info(
        "\nOpen the file and read it before you share it. The freetext_fields section "
        "contains the values themselves, and labDataName in particular is free text "
        "that could name a person. Delete anything that should not leave your site."
    )


if __name__ == "__main__":
    main()
