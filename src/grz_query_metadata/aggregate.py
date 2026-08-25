"""Aggregate per-GRZ survey reports into one spreadsheet.

Takes the grz-survey-*.json files produced by the survey command and writes an
OpenDocument spreadsheet (.ods) with one sheet per surveyed field: values down
the rows, GRZs across the columns, totals and share of all observations.

Optionally checks the free-text kit and instrument fields against a proposed
controlled vocabulary, so you can see what share of real values a candidate
enum would actually cover.

Usage
-----
    grz-aggregate-surveys grz-survey-*.json -o schema-usage.ods

    # also measure coverage of the vocabularies proposed in
    # https://github.com/jblesch/MVGenomseq/pull/1 (GRZ/vocabularies/ on its dev branch)
    grz-aggregate-surveys grz-survey-*.json -o schema-usage.ods \\
        --vocab labData.sequencerModel=GRZ/vocabularies/instrument-model.json \\
        --vocab labData.libraryPrepKit=GRZ/vocabularies/library-preparation-kit-retail-name.json
"""

from __future__ import annotations

import argparse
import collections
import json
import logging
import re
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Any

from odf.opendocument import OpenDocumentSpreadsheet

from . import __version__, ods, setup_cli_logging
from .fields import ENUM_FIELDS, FREETEXT_FIELDS
from .schema import declared_enum, enum_paths

log = logging.getLogger(__name__)


def normalise(v: str) -> str:
    """Loose key for comparing free text against a controlled value."""
    return re.sub(r"[^a-z0-9]+", "", v.lower())


def row_order(value: str) -> tuple[int, int, str]:
    """Tie-break between rows carrying the same total.

    Alphabetical, except that a wholly numeric label sorts numerically. The
    derived `duplicate_initial_submissions_per_LE` distribution is keyed by an
    exact count, and its rows tie constantly — alphabetically it would put 10
    before 2.
    """
    if re.fullmatch(r"[0-9]+", value):
        return (0, int(value), "")
    return (1, 0, value)


def sheet_name(label: str) -> str:
    # ODF forbids []:*?/\ in a sheet name. The 31-character cap is Excel's, kept
    # so that the file stays usable if someone opens it there.
    s = re.sub(r"[\[\]:*?/\\]", "", label)
    if not s:
        return "field"
    if len(s) <= 31:
        return s
    # A label with a separator in it is a metadata path: its tail is the field
    # name and its head a prefix shared with its neighbours, so the tail is what
    # distinguishes one from another. A derived metric has no such prefix, and
    # cutting its head leaves gibberish ("cate_initial_submissions_per_LE"), so
    # there the head is the end to keep.
    return s[-31:] if ("." in label or "/" in label) else s[:31]


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="grz-aggregate-surveys",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    ap.add_argument("reports", nargs="+", help="grz-survey-*.json files")
    ap.add_argument("-o", "--out", default="schema-usage.ods")
    ap.add_argument(
        "--vocab",
        action="append",
        default=[],
        help="FIELD=path/to/vocabulary.json, repeatable",
    )
    ap.add_argument(
        "--schema",
        help="grz-schema.json to read the declared enums from. Adds a zero "
        "row for every declared value that no GRZ ever used, and flags "
        "used values that the schema does not declare.",
    )
    return ap


def label_for(report: dict, path: str) -> str:
    """The column heading a report gets: its grz_id, which the survey reads out
    of the submissions themselves. The filename fallback only catches reports
    from before that was true."""
    return report.get("grz_id") or Path(path).stem


def load_reports(paths: Iterable[str]) -> list[dict]:
    """Read the survey reports, sorted by column heading. Each is tagged with a
    `_label` key holding that heading, which everything downstream reads.

    Two reports with the same heading are refused outright: they would merge
    into one column in every field sheet while the summary counted both, so
    the spreadsheet would contradict itself.
    """
    reports = []
    for p in paths:
        with open(p, encoding="utf-8") as fh:
            r = json.load(fh)
        r["_label"] = label_for(r, p)
        reports.append(r)
    reports.sort(key=lambda r: r["_label"])
    labels = [r["_label"] for r in reports]
    duplicated = sorted({label for label in labels if labels.count(label) > 1})
    if duplicated:
        sys.exit(
            f"more than one report labelled {', '.join(duplicated)} — the columns would "
            f"silently merge. Keep one report per GRZ, or rename the extra files "
            f"(the filename is the label when a report has no grz_id)."
        )
    return reports


def load_vocabs(specs: Iterable[str]) -> dict[str, set[str]]:
    """Read the FIELD=path/to/vocabulary.json specs into normalised value sets."""
    vocabs = {}
    for spec in specs:
        field, sep, path = spec.partition("=")
        if not sep or not field or not path:
            sys.exit(f"--vocab expects FIELD=path/to/vocabulary.json, got: {spec!r}")
        with open(path, encoding="utf-8") as fh:
            enum = json.load(fh).get("enum")
        if not isinstance(enum, list):
            sys.exit(f'no "enum" array in {path}')
        vocabs[field] = {normalise(v) for v in enum}
    return vocabs


def load_declared(schema_path: Path) -> dict[str, list[str]]:
    """The enum declared in the schema at each surveyed path, where there is one.

    A reference the walker cannot follow (renamed $def, missing vocabulary
    file, reference cycle) downgrades that one field to "not declared", with a
    warning — it must not abort the whole aggregation.
    """
    with open(schema_path, encoding="utf-8") as fh:
        schema_root = json.load(fh)
    declared = {}
    for label, mpath in ENUM_FIELDS.items():
        try:
            found = declared_enum(schema_root, schema_path.parent, mpath)
        except (KeyError, IndexError, ValueError, OSError, json.JSONDecodeError) as e:
            log.warning("could not walk the schema for %s: %r", label, e)
            found = None
        if found:
            declared[label] = found
    absent = sorted(set(ENUM_FIELDS) - set(declared))
    if absent:
        log.info(
            "note: no enum declared in %s for %d surveyed field(s): %s",
            schema_path.name,
            len(absent),
            ", ".join(absent),
        )

    # The reverse check: enums the schema declares at paths the survey never
    # reads. Without this, a schema that grows a new enum field would leave the
    # spreadsheet looking complete while silently omitting it.
    surveyed = set(ENUM_FIELDS.values()) | set(FREETEXT_FIELDS.values())
    try:
        unsurveyed = sorted(set(enum_paths(schema_root, schema_path.parent)) - surveyed)
    except (KeyError, IndexError, ValueError, OSError, json.JSONDecodeError) as e:
        log.warning("could not scan %s for unsurveyed enums: %r", schema_path.name, e)
        unsurveyed = []
    if unsurveyed:
        log.warning(
            "%s declares enums at %d path(s) the survey does not cover: %s",
            schema_path.name,
            len(unsurveyed),
            ", ".join(unsurveyed),
        )
    return declared


@dataclass
class Field:
    """One surveyed field, as it appears across every GRZ report.

    Everything the field's sheet and its index row need is derived here, so
    that the two cannot disagree about how many values were used or how many
    of the declared ones never occurred. The instance is never mutated after
    :func:`collect_fields` builds it, which is what makes the cached
    properties safe.
    """

    section: str  # enum_fields, freetext_fields or derived
    label: str
    per_grz: dict[str, dict[str, int]]  # GRZ label -> value -> times seen
    declared: list[str] | None = None  # enum declared in the schema, if any
    vocab: set[str] | None = None  # normalised proposed vocabulary; None = --vocab not given
    sheet: str = ""  # sheet name, assigned (deduplicated) by collect_fields

    @cached_property
    def grzs(self) -> list[str]:
        """The GRZs, in the order their columns appear."""
        return list(self.per_grz)

    @cached_property
    def merged(self) -> collections.Counter[str]:
        """All GRZs' counts folded together — the one pass everything else reads."""
        merged: collections.Counter[str] = collections.Counter()
        for counts in self.per_grz.values():
            merged.update(counts)
        return merged

    @cached_property
    def observed(self) -> set[str]:
        """Every value that actually arrived in some GRZ's submissions."""
        return set(self.merged)

    @cached_property
    def values(self) -> list[str]:
        """The rows of the sheet: every observed value plus every declared one,
        most frequent first, ties broken by :func:`row_order`. Declared values
        nobody used are included deliberately — that absence is the whole point
        of the survey, and it is invisible otherwise."""
        return sorted(
            self.observed | set(self.declared or []),
            key=lambda v: (-self.total_for(v), row_order(v)),
        )

    @cached_property
    def total(self) -> int:
        """Every observation of this field, across all GRZs."""
        return sum(self.merged.values())

    @cached_property
    def unused(self) -> list[str]:
        """Declared by the schema, never submitted by any LE: candidates for removal."""
        return [v for v in (self.declared or []) if v not in self.observed]

    @cached_property
    def covered(self) -> int:
        """Observations whose value the proposed vocabulary already contains."""
        return sum(self.total_for(v) for v in self.observed if self.in_vocab(v))

    def total_for(self, value: str) -> int:
        return self.merged[value]

    def in_vocab(self, value: str) -> bool:
        # `is not None`, never truthiness: an EMPTY proposed vocabulary is a
        # real input (it covers nothing) and must still produce the column,
        # the highlights and a 0% coverage figure.
        return self.vocab is not None and normalise(value) in self.vocab


def collect_fields(
    reports: list[dict],
    vocabs: dict[str, set[str]] | None = None,
    declared: dict[str, list[str]] | None = None,
) -> list[Field]:
    """Every field any report mentions, one sheet's worth each, ordered by
    section and then alphabetically.

    Sheet names are assigned here so that two labels which truncate to the
    same 31 characters still get distinct sheets — an ODS with two same-named
    tables is invalid, and readers would silently show only one of them.
    """
    vocabs = vocabs or {}
    declared = declared or {}
    fields = []
    taken: set[str] = set()
    for section in ("enum_fields", "freetext_fields", "derived"):
        for label in sorted({k for r in reports for k in r.get(section, {})}):
            name = sheet_name(label)
            n = 2
            while name in taken:
                suffix = f"~{n}"
                name = sheet_name(label)[: 31 - len(suffix)] + suffix
                n += 1
            taken.add(name)
            fields.append(
                Field(
                    section=section,
                    label=label,
                    per_grz={
                        r["_label"]: r.get(section, {}).get(label, {}).get("values", {}) for r in reports
                    },
                    declared=declared.get(label),
                    vocab=vocabs.get(label),
                    sheet=name,
                )
            )
    return fields


def duplicate_check_status(report: dict) -> str:
    """Whether this GRZ's database had the columns the duplicate check reads.

    A report from before the check existed carries no such key, and "unknown" is
    the honest rendering: like a "no", its duplicate figures are absent rather
    than zero.
    """
    available = report.get("duplicate_check_available")
    return "unknown" if available is None else ("yes" if available else "no")


def summary_rows(reports: list[dict]) -> list[ods.Row]:
    """The **summary** sheet: one row per GRZ report, saying how many
    submissions its database held, how many of those carried metadata and how
    many could not be parsed, plus the version that produced it.

        GRZ | submissions in table | with metadata | unparseable |
        duplicate check ran | id from | script version | generated

    "id from" says where the GRZ column heading came from: the id the
    submissions themselves record, or a `--grz-id` someone typed. Only the
    latter can have been mistyped, so it is worth seeing at a glance.

    "duplicate check ran" says whether that database had the basic_qc_passed and
    pseudonym columns the duplicate initial submission check reads. Where it did
    not, the check counted nothing, and its sheet shows that GRZ a column of
    zeroes that would otherwise read as a clean site.

    Ends with a TOTAL row, then a highlighted warning row for each of two
    things worth not missing: reports produced by different versions of this
    tool, whose columns are not strictly comparable, and GRZs whose duplicate
    check never ran.
    """
    rows = [
        ods.Row(
            [
                "GRZ",
                "submissions in table",
                "with metadata",
                "unparseable",
                "duplicate check ran",
                "id from",
                "script version",
                "generated",
            ],
            style=ods.HEAD,
        )
    ]
    rows += [
        ods.Row(
            [
                r["_label"],
                r["submissions_in_table"],
                r["submissions_with_metadata"],
                r["submissions_unparseable"],
                duplicate_check_status(r),
                r.get("grz_id_source", "unknown"),
                r["script_version"],
                r["generated"],
            ]
        )
        for r in reports
    ]
    rows.append(ods.Row())
    rows.append(
        ods.Row(
            [
                "TOTAL",
                sum(r["submissions_in_table"] for r in reports),
                sum(r["submissions_with_metadata"] for r in reports),
                sum(r["submissions_unparseable"] for r in reports),
            ],
            styles={0: ods.HEAD},
        )
    )

    versions = {r["script_version"] for r in reports}
    if len(versions) > 1:
        rows.append(ods.Row())
        rows.append(
            ods.Row(
                [f"WARNING: reports produced by different script versions: {sorted(versions)}"],
                style=ods.PROBLEM,
            )
        )

    # Without this the reader cannot tell a site with no duplicates from a site
    # where the check never ran: both show a column of zeroes.
    unchecked = [r["_label"] for r in reports if duplicate_check_status(r) != "yes"]
    if unchecked:
        rows.append(ods.Row())
        rows.append(
            ods.Row(
                [
                    f"WARNING: duplicate initial submissions were not checked for "
                    f"{', '.join(unchecked)} — that database is missing basic_qc_passed or "
                    f"pseudonym, or the report predates the check. Read their zeroes as "
                    f"'not checked', not 'none found'."
                ],
                style=ods.PROBLEM,
            )
        )
    return rows


def index_rows(fields: list[Field]) -> list[ods.Row]:
    """The **index** sheet: one row per surveyed field, as a way in to the rest
    of the file.

        section | field | values used | declared in schema |
        declared but NEVER used | sheet

    The last three columns only carry a figure where a schema was given. A
    non-zero "declared but NEVER used" count is highlighted: those are the
    fields where the schema declares vocabulary that never occurs in any GRZ's
    submissions.
    """
    rows = [
        ods.Row(
            ["section", "field", "values used", "declared in schema", "declared but NEVER used", "sheet"],
            style=ods.HEAD,
        )
    ]
    for f in fields:
        declared_count = len(f.declared) if f.declared else ""
        unused_count = len(f.unused) if f.declared else ""
        rows.append(
            ods.Row(
                [f.section, f.label, len(f.observed), declared_count, unused_count, f.sheet],
                styles={4: ods.UNUSED} if f.declared and f.unused else {},
            )
        )
    return rows


def field_header(field: Field) -> list[str]:
    """The columns of a field sheet. The last two only appear when the
    corresponding option was given, which is why the sheets are not all the
    same width."""
    header = ["value", *field.grzs, "TOTAL", "share"]
    if field.vocab is not None:
        header.append("in proposed vocabulary")
    if field.declared:
        header.append("declared in schema")
    return header


def field_rows(field: Field) -> list[ods.Row]:
    """One field's sheet: values down the rows, GRZs across the columns.

        value | <one column per GRZ> | TOTAL | share
              | in proposed vocabulary   (only with --vocab)
              | declared in schema       (only with --schema)

    A row is highlighted when it is declared by the schema but never submitted
    by any LE (its TOTAL is 0), or when it falls outside the proposed vocabulary
    (its "in proposed vocabulary" cell reads NO). The colour is emphasis only:
    both cases are readable from the columns alone.

    Below the table come a TOTAL row and, where the corresponding option was
    given, the list of declared-but-unused values and the share of observations
    the proposed vocabulary would cover.
    """
    rows = [ods.Row(field_header(field), style=ods.HEAD)]

    grzs, grand_total = field.grzs, field.total
    never_used = set(field.unused)  # membership via the Field helper, not re-derived
    for value in field.values:
        counts = [field.per_grz[g].get(value, 0) for g in grzs]
        total = field.total_for(value)
        cells: list[Any] = [
            value,
            *counts,
            total,
            round(total / grand_total, 4) if grand_total else 0,
        ]
        if field.vocab is not None:
            cells.append("yes" if field.in_vocab(value) else "NO")
        if field.declared:
            # Not an error: submissions are validated when they arrive, so a
            # value missing from today's schema simply predates a change.
            cells.append("yes" if value in field.declared else "no - predates this version")

        if value in never_used:
            style = ods.UNUSED
        elif field.vocab is not None and not field.in_vocab(value):
            style = ods.OUTSIDE_VOCAB
        else:
            style = None
        rows.append(ods.Row(cells, style=style))

    rows.append(ods.Row())
    rows.append(
        ods.Row(
            ["TOTAL", *[sum(field.per_grz[g].values()) for g in grzs], grand_total],
            styles={0: ods.HEAD},
        )
    )

    if field.declared:
        unused = field.unused
        rows.append(ods.Row())
        rows.append(
            ods.Row(
                [
                    f"declared but never used ({len(unused)} of {len(field.declared)})",
                    ", ".join(unused) if unused else "none - every declared value occurs",
                ],
                styles={0: ods.HEAD, 1: ods.UNUSED} if unused else {0: ods.HEAD},
            )
        )

    if field.vocab is not None:
        rows.append(ods.Row())
        rows.append(
            ods.Row(
                [
                    "coverage by proposed vocabulary",
                    f"{field.covered}/{field.total} observations",
                    f"{round(100 * field.covered / field.total, 1) if field.total else 0}%",
                ],
                styles={0: ods.HEAD},
            )
        )
    return rows


def build_spreadsheet(reports: list[dict], fields: list[Field]) -> OpenDocumentSpreadsheet:
    """The whole file: a summary sheet, an index sheet, then one sheet per
    surveyed field in index order."""
    doc = ods.new_document()
    ods.write_sheet(doc, "summary", summary_rows(reports), [ods.WIDE] + [ods.NARROW] * 7)
    ods.write_sheet(doc, "index", index_rows(fields), [ods.NARROW, ods.WIDE, *[ods.NARROW] * 3, ods.WIDE])
    for f in fields:
        widths = [ods.WIDE] + [ods.NARROW] * (len(field_header(f)) - 1)
        ods.write_sheet(doc, f.sheet, field_rows(f), widths)
    return doc


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    setup_cli_logging()

    reports = load_reports(args.reports)
    vocabs = load_vocabs(args.vocab)
    schema_path = Path(args.schema) if args.schema else None
    declared = load_declared(schema_path) if schema_path else {}

    fields = collect_fields(reports, vocabs, declared)

    # A --vocab label that matches no surveyed field would be silently ignored,
    # and the console would still claim its coverage was checked. Refuse instead.
    unmatched = sorted(set(vocabs) - {f.label for f in fields})
    if unmatched:
        sys.exit(
            f"--vocab field(s) not found in any report: {', '.join(unmatched)}\n"
            f"(field labels are case-sensitive; see fields.py for the surveyed labels)"
        )

    build_spreadsheet(reports, fields).save(args.out)

    log.info("wrote %s: %d field sheets from %d GRZ report(s)", args.out, len(fields), len(reports))
    for label in vocabs:
        log.info("  vocabulary coverage checked for %s", label)
    if declared and schema_path:
        unused_by_field = {f.label: f.unused for f in fields if f.declared}
        total_unused = sum(len(u) for u in unused_by_field.values())
        log.info(
            "  checked %d enum(s) against %s: %d declared value(s) never used",
            len(declared),
            schema_path.name,
            total_unused,
        )
        for label, unused in sorted(unused_by_field.items()):
            if unused:
                log.info("    %s: %s", label, ", ".join(unused))


if __name__ == "__main__":
    main()
