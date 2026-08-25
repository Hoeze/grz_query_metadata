# grz_query_metadata

Survey which metadata values are actually used across GRZ submissions, so that
proposed changes to the [GRZ metadata schema](https://github.com/BfArM-MVH/MVGenomseq_GRZ)
can be decided on evidence rather than assumption.

The proposal these questions are aimed at is
[MVGenomseq PR #1 — *Revise GRZ metadata schema for 2.0.0*](https://github.com/jblesch/MVGenomseq/pull/1),
called **the schema proposal** throughout this README.

Each GRZ runs `grz-survey-metadata` against its own submission database and
returns a small JSON report. Those reports are combined with
`grz-aggregate-surveys` into a single spreadsheet.

Typical questions this answers:

- How often is `wxs` used? It has been discussed if wxs is really needed, also because relationship to wes is unclear
- How many submissions use an ontology other than BRENDA in `tissueOntology`?
- Is `reverse` in `sequencingLayout` ever used at all?
- Do the free-text `libraryPrepKit` and `sequencerModel` values in use today map
  onto the controlled vocabularies **proposed in the schema proposal**, and what
  share of real submissions would fail if those vocabularies were adopted as
  drafted? The unmatched values are the list of entries still to be added.

## What the survey does to your database

**It reads. It never writes.** The survey issues exactly one statement:

```sql
SELECT id, submission_metadata, basic_qc_passed, pseudonym FROM submissions
```

It then inspects the metadata JSON in Python. No JSON-path SQL is used, so the
same code runs unchanged against SQLite and PostgreSQL.

`basic_qc_passed` and `pseudonym` are the only columns read outside the metadata
document, and one check uses them:
[duplicate initial submissions](#duplicate-initial-submissions). A database old
enough to predate either is detected before the query runs and read without
them — that check then reports nothing and says so, and every other count is
unaffected.

**`pseudonym` here is not a donor pseudonym.** grz-db stores the submitter's
`localCaseId` in a column named after the other field it redacts
(`_METADATA_FIELD_TO_COLUMN` maps `local_case_id` to `pseudonym`). The donor
pseudonyms live in the separate `donors` table, which this tool never reads.

### What leaves your site

The report contains **counts of values**, plus the distinct values found in
technical fields such as kit and instrument names — that is what makes the
vocabulary comparison possible. It also carries your own
`submission.genomicDataCenterId`, which is how the report labels itself; that is
your GRZ id, not anyone's patient or institution.

The following are never read and never appear in the output:

| Never touched | |
|---|---|
| `tanG` | |
| `donorPseudonym`, the `donors` table's `pseudonym` | not the `submissions.pseudonym` column, which holds `localCaseId` |
| `clinicalDataNodeId` | the KDK, not us |
| `filePath`, `fileChecksum` | |
| any date field | |
| submission ids | the database's own row ids |

Two more are read but **never written out**:

| Read, never emitted | Where from | Why it is read |
|---|---|---|
| `submitterId` | the metadata document | identifies the submitting institution under §293 SGB V |
| `localCaseId` | the `submissions.pseudonym` column | the LE's own case number |

Both are needed by one check and one check only: whether an LE sent the same
case as an `initial` submission more than once. The pair is used to group rows
in memory and is then dropped. What reaches the report is the distribution
alone — how many LEs duplicated nothing, how many duplicated once, and so on —
with no identifier and no per-institution row.

The **console** is the exception, deliberately: it names each duplicated
`submitterId` and `localCaseId`, because it runs at the site that can go and
look those cases up. That output is local. Do not paste it into a ticket or a
mail to the coordination — send the report file, which cannot carry them. See
[duplicate initial submissions](#duplicate-initial-submissions).

**Open the report and read it before you send it.** The file is small,
human-readable JSON, and the only part that is not a plain count is the
`freetext_fields` section. One field there deserves particular attention:
`labDataName` is free text and could in principle contain something identifying,
for example `"Blut DNA Müller"`. Delete anything that should not leave your site
— removing values from the JSON does not break the aggregation. The tool prints
this reminder on every run.

## Installation

This is a Python package managed with [uv](https://docs.astral.sh/uv/). uv is
the only thing you need to install yourself; it fetches a suitable Python
(3.12+) and the dependencies (`sqlalchemy`, `psycopg[binary]`, `PyYAML`,
`odfpy`) into an isolated environment, so nothing is added to your `grz_tools`
installation. The PostgreSQL driver is included deliberately — SQLite works out
of the box either way, but a site whose submission database is PostgreSQL must
not discover a missing driver at the moment it connects — and it is psycopg 3
with the same version pin as `grz-db`, so the sites see one driver ecosystem.
Plain `postgresql://` URLs are routed to it automatically.

The quickest way to run a pinned version without installing anything
permanently:

```bash
uvx --from git+https://github.com/jblesch/grz_query_metadata@v2.0.0 grz-survey-metadata --help  # x-release-please-version
```

To install the two commands onto your PATH instead:

```bash
uv tool install git+https://github.com/jblesch/grz_query_metadata@v2.0.0  # x-release-please-version
```

(The `x-release-please-version` markers keep the pinned tag current: every
release bumps these lines automatically.)

Or from a clone, which is also what you want when working on the code:

```bash
git clone https://github.com/jblesch/grz_query_metadata
cd grz_query_metadata
uv sync
uv run grz-survey-metadata --help
```

Both commands print the package version with `--version`, and that same version
is recorded in every report, so a report can always be traced back to the tag
that produced it.

## Running the survey

Point it at your database. That is the whole of it:

```bash
# SQLite: a plain path is fine
grz-survey-metadata --db-url /path/to/submission.db.sqlite

# PostgreSQL
grz-survey-metadata --db-url postgresql://user@host/grzdb

# or take the URL from your grz config file
grz-survey-metadata --config-file /etc/grz/config.yaml
```

Prefix these with `uv run` (in a clone) or `uvx --from git+https://github.com/jblesch/grz_query_metadata`
(without installing) if you did not `uv tool install`.

This writes `grz-survey-GRZK00123-<date>.json` in the working directory.

| Flag | Effect |
|---|---|
| `--db-url URL` | the database, as a SQLAlchemy URL or a path to a SQLite file |
| `--config-file FILE` | take the URL from a grz config file instead |
| `--out FILE` | write somewhere other than the default filename |
| `--grz-id ID` | override the site id — see below, normally unnecessary |

### Where the site id comes from

You do not type it. Every submission records the GRZ that received it in
`submission.genomicDataCenterId` (`GRZXXXnnn`, e.g. `GRZK00123`), so the survey
reads it out of the data, names the report after it and labels your column with
it. A retyped id is the one part of the report that could be wrong; a derived
one cannot be.

The run prints what it concluded, and the report records it under
`grz_id_source`, which the aggregated spreadsheet shows in an **id from** column
so the coordinator can see which reports were labelled by hand:

```
reporting as GRZK00123 (from submission/genomicDataCenterId)
```

`--grz-id` exists for the two cases where the data cannot answer. Both stop the
run rather than guessing:

| Situation | What happens |
|---|---|
| no submission records an id | error; pass `--grz-id` to label the report yourself |
| submissions from **more than one** GRZ | error listing each id and how many submissions carry it; pass `--grz-id` to say which this report is for |

The second is worth a look before you reach for the flag: a handful of
submissions bearing another GRZ's id in your database is a finding in itself,
which is why the error prints the counts rather than silently picking the
majority.

Every distinct value of every surveyed field is reported, with no cap. A
free-text field such as `labDataName` can hold thousands of one-off values, so
the report is not always small — but a truncated one would understate exactly
the long tail a controlled vocabulary has to account for.


## Aggregating the reports

Run by whoever collects the reports:

```bash
grz-aggregate-surveys grz-survey-*.json -o schema-usage.ods
```

The output is an OpenDocument spreadsheet, which LibreOffice, Excel, Numbers and
Google Sheets all open, and which `pandas.read_excel(..., engine="odf")` reads
directly. Counts are written as numbers rather than text, so sorting and summing
in the sheet work as expected.

To also measure how well a proposed controlled vocabulary would cover the
free-text values actually in use:

```bash
grz-aggregate-surveys grz-survey-*.json -o schema-usage.ods \
    --vocab labData.sequencerModel=path/to/instrument-model.json \
    --vocab labData.libraryPrepKit=path/to/library-preparation-kit-retail-name.json
```

Both vocabulary files come from the schema proposal itself — they are added by
[PR #1](https://github.com/jblesch/MVGenomseq/pull/1) and live in
`GRZ/vocabularies/` on its `dev` branch:

| File | Compared against |
|---|---|
| [`GRZ/vocabularies/instrument-model.json`](https://github.com/jblesch/MVGenomseq/blob/dev/GRZ/vocabularies/instrument-model.json) | the free-text `sequencerModel` values in use today |
| [`GRZ/vocabularies/library-preparation-kit-retail-name.json`](https://github.com/jblesch/MVGenomseq/blob/dev/GRZ/vocabularies/library-preparation-kit-retail-name.json) | the free-text `libraryPrepKit` values in use today |

So the coverage figure answers a question about the proposal specifically: if
these two fields became closed enums as drafted, what share of what the
Leistungserbringer have already submitted to the GRZs would still validate.

### Finding the values nobody uses

The survey counts what it sees, so an enum value that never occurs is simply
absent from the report. That is a problem, because "never used by anyone" is
precisely the finding that justifies dropping a value. Pass the schema and those
gaps are filled in:

```bash
grz-aggregate-surveys grz-survey-*.json -o schema-usage.ods \
    --schema path/to/grz-schema.json
```

This reads the enum declared at each surveyed path — following `$ref` into
`GRZ/vocabularies/` where needed — and adds a row at zero for every declared
value no GRZ used. Each sheet then shows the complete vocabulary instead of only
the part in use. Rows are colour-coded:

| Colour | Meaning | Also readable from |
|---|---|---|
| blue | declared in the schema, **never used by anyone** — a candidate for removal | `TOTAL` is 0 |
| amber | outside the proposed vocabulary, where `--vocab` was also given | `in proposed vocabulary` is `NO` |

The colour is emphasis, never the only carrier of a finding: every highlighted
row can be found from a column instead, so the sheet keeps its meaning if it is
exported to CSV or read by a script.

A value present in the data but absent from the schema is **not** flagged as a
problem. Submissions are validated by grz-tools when they arrive, so such a value
was valid against the schema version in force at the time; the database simply
spans several versions. The `declared in schema` column records it as
`no - predates this version` and nothing is highlighted. It is still worth a
glance for one reason: if a whole file is full of them, the wrong schema file
was passed.

Every sheet ends with a `declared but never used (8 of 10)` line naming them, the
**index** sheet carries a per-field count, and the same summary goes to the
console:

```
checked 20 enum(s) against grz-schema.json: 75 declared value(s) never used
  labData.libraryType: panel, panel_lr, wes_lr, wgs, wgs_lr, wxs_lr, other, unknown
  labData.sequencingLayout: single-end, reverse, other
```

Pass the schema the submissions were actually **validated against** — that is
`GRZ/grz-schema.json` on `main`, not the schema proposal. The job here is to
measure the vocabulary in force; measuring against the proposal is what `--vocab`
does. `--schema` and `--vocab` can be given together and answer different halves
of the question: which existing values are dead, and which proposed values are
missing.

The `.ods` contains:

- **summary** — submissions seen per GRZ, where each GRZ's id came from, whether
  that database recorded QC status, and a warning row if script versions differ
  or if some GRZ never ran the duplicate check
- **index** — every surveyed field, with how many distinct values were used and —
  with `--schema` — how many values the schema declares and how many of those
  were never used
- **one sheet per field** — values down the rows, GRZs across the columns, plus
  totals and share of all observations

Where `--vocab` was given, each value is marked as in or out of the proposed
vocabulary, unmatched rows are highlighted, and the sheet ends with an overall
coverage figure such as `7/12 observations (58.3%)`. That figure is the point of
the exercise: it says whether a candidate enum is usable as drafted, and the
highlighted rows are the list of values still to be mapped or added.

## What is surveyed

**Controlled vocabularies today** — `submissionType`, `coverageType`,
`diseaseType`, `genomicStudyType`, `genomicStudySubtype`, `gender`, `relation`,
`noScopeJustification`, `sampleConservation`, `sequenceType`, `sequenceSubtype`,
`fragmentationMethod`, `libraryType`, `enrichmentKitManufacturer`,
`sequencingLayout`, `tumorCellCount.method`, `referenceGenome`, `fileType`,
`readOrder`, `checksumType`.

**Free text today** — `tissueOntology.name` and `.version`, `tissueTypeId`,
`tissueTypeName`, `libraryPrepKit`, `libraryPrepKitManufacturer`,
`sequencerModel`, `sequencerManufacturer`, `kitName`, `kitManufacturer`,
`enrichmentKitDescription`, `labDataName`, `bioinformaticsPipelineName`,
`callerUsed.name`.

**Derived checks** that a plain value count cannot answer:

| Check | Question it answers |
|---|---|
| `tissueTypeId_is_BTO_format` | how many identifiers match `^BTO:[0-9]{7}$` |
| `duplicate_initial_submissions_per_LE` | how many Leistungserbringer re-sent a case as a QC-passed `initial`, and exactly how many times each |
| `initial_submissions_checked_for_duplicates` | how much of the input the check above could use, and why the rest was left out |

### Duplicate initial submissions

A case is meant to arrive once as `initial`. Anything after that is a
`followup`, an `addition` or a `correction`. A second `initial` for a case an LE
has already sent is therefore a re-submission nobody intended, and it inflates
every count in this report along with everything else counted per submission.

**Only submissions that passed basic QC count.** One that failed was rejected,
and the LE was meant to send it again — counting that replacement as a duplicate
would report the process working as a fault. This is grz-db's own definition:
its `_qc_passed_initial_of` and its one-initial-per-case index both mean
`basic_qc_passed IS TRUE`, and `detailed_qc_passed` is deliberately not used,
since in-depth QC runs on a selected sample and filtering on it would discard
most submissions rather than the failed ones. A `NULL` — QC not yet decided —
is not a pass.

Because grz-db's index already forbids a second QC-passed initial per case,
what this check surfaces is the part that index does not reach: submissions
from before it existed, and those with no case linked.

The `localCaseId` half of the key comes from the `submissions.pseudonym` column,
not from the metadata document, because grz-db **redacts** `localCaseId` in the
JSON it stores, replacing it with `REDACTED_LOCAL_CASE_ID`. Reading the document
would key every redacted submission of an LE on that one placeholder and report
them all as duplicates of each other. A submission whose column still holds a
placeholder, or nothing, identifies no case and is counted apart rather than
grouped, which is also how grz-db's own `SubmitterLocalCaseResolver` treats it.

The check groups the QC-passed initial submissions by `(submitterId, localCaseId)`,
counts how many each LE sent beyond the first, and reports the **distribution
over LEs**:

```json
"duplicate_initial_submissions_per_LE": {
  "_total": 14, "_distinct": 4,
  "values": { "0": 11, "1": 1, "2": 1, "7": 1 }
}
```

Read as: of the 14 LEs that sent anything, 11 duplicated nothing, one
duplicated once, one twice, and one sent seven initial submissions more than it
should have. The key is the exact count, never a range and never the LE. An LE
with several duplicated cases is counted across all of them, so `7` may be one
case sent eight times or seven cases sent twice.

Neither identifier appears in the report. Note that an exact count is a sharper
key than a range: an LE with an unusual number is distinguishable by its row
alone to anyone who already knows roughly how much it submits.

`localCaseId` is the LE's own numbering, so it is only meaningful within one LE:
two LEs both using `case-1` is not a duplicate, and the check does not treat it
as one.

A companion counter accounts for every initial submission, whether it was
compared or not:

```json
"initial_submissions_checked_for_duplicates": {
  "values": {
    "yes": 812,
    "no - basic QC not passed": 37,
    "no - localCaseId redacted": 9,
    "no - no submitterId or localCaseId": 4
  }
}
```

Each reason is written even at zero, so the duplicate count can be read as a
total rather than a floor. `no - basic QC not passed` also covers submissions
whose QC is still undecided; if it swallows nearly everything, check whether
this database records QC at all.

Across GRZs, two things about this metric are worth keeping in mind. A GRZ whose
database lacks either column counts nothing, so its column is all
zeroes and reads like a clean site: the report records
`duplicate_check_available`, the **summary** sheet has a `duplicate check ran`
column, and a highlighted warning row names any GRZ whose check never ran. And the TOTAL row adds LEs
across GRZs, so an LE that submits to more than one GRZ is counted once per GRZ.
Deduplicating would need the identifiers, which deliberately never leave the
sites.

### Acting on it

The report says how many; the console says which. On every run the survey names
every duplicated case, worst first, one line per affected LE:

```
warning: 2 of 14 Leistungserbringer sent the same case as a QC-passed initial submission more than once:
  (these submitterIds and localCaseIds stay on this machine — the report carries only how many duplicates each LE had, and no identifier)
    260000001: ABC-17 (3x), ABC-92 (2x)
    260000002: XY-4 (2x)
```

`ABC-17 (3x)` means that case arrived as a QC-passed `initial` submission three
times, so two of the three were unintended. Nothing is truncated, and nothing here
goes into the report file. When no LE duplicated anything the survey says so
instead, with the number of LEs it checked.

## Adding a field

Both dictionaries in [`src/grz_query_metadata/fields.py`](src/grz_query_metadata/fields.py)
map a label to a path through the metadata document, where `[]` descends into
an array:

```python
"labData.sequencingLayout": "donors[]/labData[]/sequencingLayout",
```

Add an entry and it is counted — by the survey and, when `--schema` is given,
by the schema walk as well, since both read the same dictionary. There is no
version constant to bump: the report records the package version, which
release-please derives from the commit history (see below).

## Development

```
src/grz_query_metadata/
  fields.py     the metadata paths that are surveyed — the one place a field is declared
  survey.py     counts values in a submission database        → grz-survey-metadata
  schema.py     reads the enums declared in a GRZ JSON Schema
  aggregate.py  builds the spreadsheet from the JSON reports  → grz-aggregate-surveys
  ods.py        a thin writer over odfpy: sheets, rows, cell styles
tests/
```

Inside `aggregate.py`, each sheet is produced by one function whose docstring
states exactly what that sheet contains and which columns appear under which
option:

| Function | Produces |
|---|---|
| `summary_rows(reports)` | the **summary** sheet |
| `index_rows(fields)` | the **index** sheet |
| `field_rows(field)` | one field's sheet, `field_header(field)` being its columns |
| `build_spreadsheet(reports, fields)` | all of the above, in order |

They return lists of `ods.Row`, not spreadsheet objects, so what a sheet
contains can be read — and tested — without a spreadsheet in the picture. The
arithmetic they share lives on `Field`, so a sheet and its index row cannot
disagree about how many values were used or how many declared ones never
occurred.

```bash
uv sync --all-groups     # set up the environment
uv run pytest            # tests
uv run tox               # everything CI runs: format, lint, mypy, tests
uv run ruff format .     # apply formatting
```

CI runs the same tox environments on every push and pull request.

### Releasing

Versioning is handled by
[release-please](https://github.com/googleapis/release-please) from
[Conventional Commits](https://www.conventionalcommits.org/), so commit messages
decide the next version:

| Commit prefix | Effect |
|---|---|
| `fix: ...` | patch release |
| `feat: ...` | minor release |
| `feat!: ...` or a `BREAKING CHANGE:` footer | major release |
| `chore: ...`, `docs: ...`, `test: ...` | no release |

Merging to `main` opens or updates a release PR that bumps the version in
`pyproject.toml` and writes `CHANGELOG.md`. Merging that PR tags the release
(`v1.2.0`), which is what sites should install and what the `script_version`
field of every report then reports back. Nothing is published to PyPI; installs
go straight from the git tag.
