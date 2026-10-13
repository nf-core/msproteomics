#!/usr/bin/env python3
"""File-by-file comparator for golden FragPipe differential tests.

Implements exactly the normalization rules measured in FragPipe determinism runs
(FragPipe 24.0 run twice on the same input) and nothing else. Every other byte must be identical.
Rules are chosen by the file's path relative to the compared root (never by
content). There are two rule profiles because the measurements differ:

``dda`` (default; TMT label check, DDA LFQ, glyco DDA; doc sections "Result"
and "Normalization rules for golden tests"):

=======================================  ======================================
File                                     Normalization before byte comparison
=======================================  ======================================
``log_*.txt`` (FragPipe log), ``?/``     SKIP
``*.pepindex``                           drop the first 8 bytes, path prefixes
                                         -> ``<RUN>``
``*.prot.xml``                           tmp-UUID dir -> ``<TMP>/``, path
                                         prefixes -> ``<RUN>``, ISO
                                         ``date=``/``time=``
``*.pepXML``, ``*.pep.xml``              path prefixes, ISO ``date=``/``time=``,
                                         tmp-UUID dir -> ``<TMP>/`` (only where the
                                         caller gives tmp roots; iProphet's
                                         ``combined.pep.xml`` carries the
                                         ``philosopher --temp`` dir in
                                         ``summary_xml``)
``filter.log`` (Philosopher)             paths, ``time=`` stripped, Parsing lines
                                         as a set, counts within the measured
                                         envelope (``philosopher_tables``)
``psm.tsv``                              Philosopher envelope: Qvalue bound,
                                         <=1 boundary row (``philosopher_tables``)
``peptide/ion/protein.tsv``              spectral counts +- boundary PSMs, and the
                                         IonQuant feature columns (Intensity,
                                         Retention, Ion Mobility) of an entry a
                                         boundary PSM maps to;
                                         peptide.tsv also Prev/Next AA dropped and
                                         the Qvalue rule (measured Philosopher
                                         variation: Next AA 5/6 vs 1/6 runs)
``combined_{protein,peptide,``           byte for byte, except rows of entries a
``modified_peptide,ion}.tsv``,           boundary PSM (sibling <group>/psm.tsv,
``msstats{,_ptm}.csv``,                  within the psm.tsv envelope) maps to: those
``reprint.{spc,int}.tsv``                may differ in that group's spectral-count,
                                         intensity, RT FWHM and ion-mobility columns
``FPOP_peptides.tsv``,                   byte for byte except the Peptide Sequence
``FPOP_sites.tsv``                       column, compared as a set per row (python set join)
``*.workflow``                           path prefixes -> ``<RUN>``, the ``# DIA-NN
                                         version`` line case-folded (Beta/beta: host
                                         FragPipe config, TabConfig.java:566-584)
``fragger*.params``,                     path prefixes -> ``<RUN>``
``*.job``, ``*.fp-manifest``,
``filelist_*.txt``,
``experiment_annotation.tsv``,
``DbOptimizerConfig.yaml`` (FP-Meta: its projectDir/fastaFile/outdir are the run's paths)
everything else                          raw bytes
=======================================  ======================================

peptide.tsv/ion.tsv get the full envelope whenever Philosopher produced them and no
IonQuant step rewrote them (decided per file from the manifest by ``golden_step_check``;
without a manifest the profile decides: ``dia`` full, ``dda`` counts).

``dia`` (doc sections "DIA section" and "DIA section 2") = ``dda`` plus:

=======================================  ======================================
``MSBooster/spectraRT.predicted.bin``    SKIP
``*.pdf``                                ``/CreationDate (D:...)`` -> constant
``report.log.txt`` (DIA-NN)              path prefixes, drop the
                                         ``Current date and time:`` text,
                                         ``[m:ss] `` prefixes and the host's
                                         ``Logical CPU cores:`` count
``MSBooster/spectraRT*.tsv``             sort rows (header first)
``peptide.tsv``, ``ion.tsv``             full Philosopher envelope (Qvalue
                                         rule; ion.tsv keyed on stripped
                                         sequence+Charge+M/Z, Expectation/Prev
                                         AA/Next AA dropped; peptide Prev/Next AA)
``msbooster_params.txt``                 path prefixes -> ``<RUN>``
``easypqp_files/cmds.txt``               path prefixes -> ``<RUN>``, then every
                                         ``<RUN>/<dirs>/<name>`` -> ``<RUN>/<name>``
                                         (gen_con_spec_lib.py resolves every input
                                         path, :311-312, :451, :555, :566-568, so a
                                         staged symlink names its producer's dir)
``dia-quant-output/report.tsv``,         path prefixes -> ``<RUN>`` (DIA-NN writes
``report.{pg,pr,gg,unique_genes}``       each run's ``--f`` path, CmdDiann.java:389-397,
``_matrix.tsv``, ``report.stats.tsv``    as File.Name / run column)
``dia-quant-output/abundance_*.tsv``     path prefixes -> ``<RUN>``, then the run
                                         columns (header ``<RUN>/...``, contiguous at
                                         the end) sorted by name in every row:
                                         FragReporter 1.1.2 orders them by a HashMap
                                         keyed by the absolute run path
                                         (DiannPrReport.runInfo/getRuns), so the order
                                         follows String.hashCode of the work dir
``*_diatracer.mzML``                     path prefixes in ``location="..."``, ``startTimeStamp``
                                         (wall clock of the run), and the index
                                         byte offsets (``<offset>``, ``<indexListOffset>``)
                                         that shift with the header length; every
                                         spectrum line is compared exactly
``*_sub.mzML``                           as ``*_diatracer.mzML``: WriteSubMzml's MZMLWriter
                                         writes the LC-MS file's directory as the
                                         ``sourceFile`` location (WriteSubMzml.java:217-233),
                                         so the index offsets shift with its length
=======================================  ======================================

The path rewrite inside ``filter.log`` and ``.pepindex`` was measured: the
determinism scripts (``cmpdirs.py``, ``cmp_runs.py``) rewrite the run dir
before every other rule, and ``.pepindex`` headers embed the FASTA path (a
different path length also changes the 4-byte header length).

Path prefixes are given per side (``--expected-prefix`` / ``--actual-prefix``)
and replaced longest first, only at a path boundary. Missing and extra files
are failures, and so is comparing nothing (no file other than SKIPs).
Exit code: 0 all MATCH/SKIP, 1 any DIFF/MISSING/EXTRA/EMPTY.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path, PurePosixPath

sys.path.insert(0, str(Path(__file__).resolve().parent))
import philosopher_tables as pt  # noqa: E402

RUN_TOKEN = b"<RUN>"
TMP_TOKEN = b"<TMP>/"
TIMESTAMP_TOKEN = b"<TIMESTAMP>"
PEPINDEX_HEADER_BYTES = 8  # 4 header bytes + 4-byte header length (determinism doc)
CHUNK = 1 << 20
PROFILES = ("dda", "dia")

ISO_ATTR_RE = re.compile(rb'\b(date|time)="\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}"')
UUID_DIR = rb"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}/"
MZML_START_RE = re.compile(rb'startTimeStamp="[^"]*"')
MZML_OFFSET_RE = re.compile(rb"(<offset idRef=\"[^\"]*\">|<indexListOffset>)\d+(?=</)")
PDF_CREATION_RE = re.compile(rb"/CreationDate *\(D:[0-9+'Z-]+\)")
DIANN_DATE_RE = re.compile(rb"(?m)^Current date and time: .*$")
DIANN_ELAPSED_RE = re.compile(rb"(?m)^\[\d+:\d{2}\] ")
# DIA-NN logs the host's core count (PXD047839: 128 on the golden host, 96 on Tower), not a run setting
DIANN_CORES_RE = re.compile(rb"(?m)^Logical CPU cores: \d+$")

SKIP = "skip"
PEPINDEX = "path+skip-first-8-bytes"
PROTXML = "path+tmp-uuid+iso-timestamp"
PEPXML = "path+iso-timestamp"
PATH_ONLY = "path"
# FragPipe's saved workflow: "# DIA-NN version 1.8.2 Beta 8" (Diann.fallBackDiannVersion, no DIA-NN binary
# configured) or "beta 8" (the binary's --version, binary configured in the workflow or the host's FragPipe
# cache; TabConfig.java:566-584): host state, not a run setting
WORKFLOW_FILE = "path+diann-version-case"
DIANN_VERSION_LINE_RE = re.compile(rb"(?m)^# DIA-NN version .*$")
RAW = "raw"
# Philosopher filter/report nondeterminism (determinism doc section of that name)
PHILO_FILTER_LOG = "philosopher-filter-log-envelope"
PHILO_PSM = "philosopher-psm-envelope"
PHILO_COUNTS = "philosopher-spectral-count-envelope"
PHILO_FULL = "philosopher-peptide-ion-envelope"
PHILO_RUN = "philosopher-boundary-run-table"
# FPOP Quant tables: the peptide list of a group is the join of a Python set
FPOP_SET_COLUMN = "fpop-peptide-list-as-set"
# DIA profile only
PDF = "pdf-creationdate"
DIANN_LOG = "path+diann-date+elapsed"
SORT_ROWS = "sort-rows"
EASYPQP_CMDS = "path+basename-under-run"
FRAGREPORTER_TABLE = "path+sort-run-columns"
MZML_INDEXED = "path+start-time+index-offsets"

FPOP_TABLES = ("FPOP_peptides.tsv", "FPOP_sites.tsv")
FPOP_LIST_COLUMN = "Peptide Sequence"

DIANN_TABLES = re.compile(r"^report(\.(pg|pr|gg|unique_genes)_matrix|\.stats)?\.tsv$")
FRAGREPORTER_TABLES = re.compile(r"^abundance_.*\.tsv$")
RUN_SUBDIRS_RE = re.compile(re.escape(RUN_TOKEN) + rb"/(?:[^/\s'\"]+/)+")
PATH_ONLY_NAMES = re.compile(
    r"^(fragger(_[a-z_]+)?\.params|.*\.workflow|.*\.job|.*\.fp-manifest|filelist_.*\.txt|msbooster_params\.txt|experiment_annotation\.tsv|DbOptimizerConfig\.yaml)$"
)


def rule_for(rel: str, profile: str = "dda", ionquant_rewritten: bool | None = None) -> str:
    """Normalization rule for a file, from its path relative to the compared root.

    ``ionquant_rewritten`` says whether an IonQuant step rewrote peptide.tsv/ion.tsv
    (known from the manifest: a later step outputs the same path). Philosopher-produced
    and not rewritten -> full envelope; rewritten -> spectral counts only. ``None``
    (no manifest, e.g. the ``dirs``/``pairs`` CLI) falls back to the profile:
    ``dia`` runs no IonQuant, ``dda`` does.
    """
    if profile not in PROFILES:
        raise ValueError(f"unknown rule profile {profile!r}; use one of {PROFILES}")
    parts = PurePosixPath(rel).parts
    name = parts[-1]
    if "?" in parts[:-1] or re.fullmatch(r"log_.*\.txt", name):
        return SKIP
    if profile == "dia":
        in_msbooster = "MSBooster" in parts[:-1]
        if in_msbooster and name == "spectraRT.predicted.bin":
            return SKIP
        if name.endswith(".pdf"):
            return PDF
        if name == "report.log.txt":
            return DIANN_LOG
        if in_msbooster and re.fullmatch(r"spectraRT.*\.tsv", name):
            return SORT_ROWS
        if name == "msbooster_params.txt":
            return PATH_ONLY
        if name == "cmds.txt" and "easypqp_files" in parts[:-1]:
            return EASYPQP_CMDS
        if "dia-quant-output" in parts[:-1] and DIANN_TABLES.match(name):
            return PATH_ONLY
        if "dia-quant-output" in parts[:-1] and FRAGREPORTER_TABLES.match(name):
            return FRAGREPORTER_TABLE
    # Philosopher-produced tables and log, both profiles (doc section "Philosopher
    # filter/report nondeterminism"; "Workflow scope": TMT/DDA LFQ get the psm.tsv
    # and spectral-count rules because IonQuant rewrites peptide.tsv/ion.tsv).
    if name == "filter.log":
        return PHILO_FILTER_LOG
    if name == "psm.tsv":
        return PHILO_PSM
    if name in ("peptide.tsv", "ion.tsv"):
        rewritten = profile == "dda" if ionquant_rewritten is None else ionquant_rewritten
        return PHILO_COUNTS if rewritten else PHILO_FULL  # doc "Workflow scope"
    if name == "protein.tsv":
        return PHILO_COUNTS
    if name in pt.RUN_TABLES and "dia-quant-output" not in parts[:-1]:
        return PHILO_RUN
    if name in FPOP_TABLES:
        return FPOP_SET_COLUMN
    if name.endswith(".pepindex"):
        return PEPINDEX
    if name.endswith(".prot.xml"):
        return PROTXML
    if name.endswith((".pepXML", ".pep.xml")):
        return PEPXML
    if name.endswith(".workflow"):
        return WORKFLOW_FILE
    if PATH_ONLY_NAMES.match(name):
        return PATH_ONLY
    if name.endswith(("_diatracer.mzML", "_sub.mzML")):
        return MZML_INDEXED
    return RAW


@dataclass
class Side:
    """Run-specific strings of one side: path prefixes and java.io.tmpdir roots."""

    prefixes: list[str] = field(default_factory=list)
    tmp_roots: list[str] = field(default_factory=list)

    def rewrite_paths(self, data: bytes) -> bytes:
        for p in sorted({x.rstrip("/") for x in self.prefixes if x}, key=len, reverse=True):
            data = re.sub(re.escape(p.encode()) + rb'(?=[/"\s<:]|$)', RUN_TOKEN, data)
        return data

    def rewrite_tmp(self, data: bytes) -> bytes:
        for t in sorted({x.rstrip("/") for x in self.tmp_roots if x}, key=len, reverse=True):
            data = re.sub(re.escape(t.encode()) + rb"/" + UUID_DIR, TMP_TOKEN, data)
        return data


def _iso(data: bytes) -> bytes:
    return ISO_ATTR_RE.sub(lambda m: m.group(1) + b'="' + TIMESTAMP_TOKEN + b'"', data)


def _split_rows(data: bytes) -> tuple[list[bytes], bool]:
    trailing = data.endswith(b"\n")
    rows = data.split(b"\n")
    if trailing:
        rows = rows[:-1]
    return rows, trailing


def _join_rows(rows: list[bytes], trailing: bool) -> bytes:
    return b"\n".join(rows) + (b"\n" if trailing else b"")


def _sort_rows(data: bytes, header: bool) -> bytes:
    rows, trailing = _split_rows(data)
    if header:
        return _join_rows(rows[:1] + sorted(rows[1:]), trailing)
    return _join_rows(sorted(rows), trailing)


def _sort_run_columns(data: bytes) -> bytes:
    """Sort the ``<RUN>/...`` columns (contiguous at the end of the header) by name in every row."""
    rows, trailing = _split_rows(data)
    header = rows[0].split(b"\t")
    runs = [i for i, cell in enumerate(header) if cell.startswith(RUN_TOKEN + b"/")]
    if not runs:
        raise ValueError("FragReporter table has no <RUN>/ run columns (path prefixes not given?)")
    first = runs[0]
    if runs != list(range(first, len(header))):
        raise ValueError(f"FragReporter run columns are not contiguous at the end of the header: {runs}")
    order = list(range(first)) + sorted(runs, key=lambda i: header[i])
    out = []
    for row in rows:
        cells = row.split(b"\t")
        if len(cells) != len(header):
            raise ValueError(f"FragReporter row has {len(cells)} cells, header {len(header)}")
        out.append(b"\t".join(cells[i] for i in order))
    return _join_rows(out, trailing)


def normalize(data: bytes, rule: str, side: Side) -> bytes:
    if rule == PROTXML:
        return _iso(side.rewrite_paths(side.rewrite_tmp(data)))
    if rule == PEPXML:
        return _iso(side.rewrite_paths(side.rewrite_tmp(data)))
    if rule == WORKFLOW_FILE:
        return DIANN_VERSION_LINE_RE.sub(lambda m: m.group(0).lower(), side.rewrite_paths(data))
    if rule == PATH_ONLY:
        return side.rewrite_paths(data)
    if rule == PEPINDEX:
        return side.rewrite_paths(data[PEPINDEX_HEADER_BYTES:])
    if rule == PDF:
        return PDF_CREATION_RE.sub(b"/CreationDate " + TIMESTAMP_TOKEN, data)
    if rule == DIANN_LOG:
        return DIANN_CORES_RE.sub(b"Logical CPU cores: <N>", DIANN_ELAPSED_RE.sub(b"", DIANN_DATE_RE.sub(b"", side.rewrite_paths(data))))
    if rule == SORT_ROWS:
        return _sort_rows(data, header=True)
    if rule == EASYPQP_CMDS:
        return RUN_SUBDIRS_RE.sub(RUN_TOKEN + b"/", side.rewrite_paths(data))
    if rule == FRAGREPORTER_TABLE:
        return _sort_run_columns(side.rewrite_paths(data))
    raise ValueError(f"rule {rule} is not a text normalization")


PHILO_RULES = {
    PHILO_FILTER_LOG: lambda pe, pa, e, a, n: pt.compare_filter_log(pe, pa, e, a),
    PHILO_PSM: lambda pe, pa, e, a, n: pt.compare_psm(pe, pa, e, a),
    PHILO_COUNTS: lambda pe, pa, e, a, n: pt.compare_entry_table(pe, pa, e, a, n, full=False),
    PHILO_FULL: lambda pe, pa, e, a, n: pt.compare_entry_table(pe, pa, e, a, n, full=True),
    PHILO_RUN: lambda pe, pa, e, a, n: pt.compare_run_table(pe, pa, e, a, n),
}


@dataclass
class Result:
    path: str
    status: str  # MATCH, DIFF, MISSING, EXTRA, SKIP, EMPTY
    rule: str
    detail: str = ""


FAIL_STATUSES = ("DIFF", "MISSING", "EXTRA", "EMPTY")


def _first_diff_bytes(a: bytes, b: bytes) -> str:
    n = min(len(a), len(b))
    off = next((i for i in range(n) if a[i] != b[i]), n)
    line = a.count(b"\n", 0, off) + 1
    la = a.split(b"\n")[line - 1] if line - 1 < a.count(b"\n") + 1 else b""
    lb = b.split(b"\n")[line - 1] if line - 1 < b.count(b"\n") + 1 else b""
    return (
        f"first difference at byte {off}, line {line} (sizes {len(a)} vs {len(b)}): "
        f"expected {la[:200]!r} actual {lb[:200]!r}"
    )


def _compare_raw_stream(pe: Path, pa: Path) -> str | None:
    """Stream-compare two files; return None or a difference description."""
    se, sa = pe.stat().st_size, pa.stat().st_size
    off = 0
    line = 1
    with pe.open("rb") as fe, pa.open("rb") as fa:
        while True:
            be, ba = fe.read(CHUNK), fa.read(CHUNK)
            if be != ba:
                n = min(len(be), len(ba))
                i = next((k for k in range(n) if be[k] != ba[k]), n)
                line += be.count(b"\n", 0, i)
                return f"first difference at byte {off + i}, line {line} (sizes {se} vs {sa})"
            if not be:
                return None
            line += be.count(b"\n")
            off += len(be)


def _normalize_mzml_line(line: bytes, side: Side) -> bytes:
    line = MZML_OFFSET_RE.sub(rb"\1<OFFSET>", MZML_START_RE.sub(b'startTimeStamp="' + TIMESTAMP_TOKEN + b'"', line))
    return side.rewrite_paths(line)


def _compare_mzml_stream(pe: Path, pa: Path, exp: Side, act: Side) -> str | None:
    """Line-wise compare of an indexed mzML (1-2 GB: no read_bytes); None or the first differing line."""
    with pe.open("rb") as fe, pa.open("rb") as fa:
        for n, (le, la) in enumerate(zip(fe, fa), 1):
            ne, na = _normalize_mzml_line(le, exp), _normalize_mzml_line(la, act)
            if ne != na:
                return f"first difference at line {n}: expected {ne[:200]!r} actual {na[:200]!r}"
        if fe.read(1) or fa.read(1):
            return "files have a different number of lines"
    return None


def compare_files(pe: Path, pa: Path, rel: str, exp: Side, act: Side, profile: str = "dda",
                  ionquant_rewritten: bool | None = None) -> Result:
    """Compare one file pair; ``rel`` (path relative to the compared root) selects the rule."""
    rule = rule_for(rel, profile, ionquant_rewritten)
    if rule == SKIP:
        return Result(rel, "SKIP", rule, "excluded by determinism rule")
    if rule == RAW:
        d = _compare_raw_stream(pe, pa)
        return Result(rel, "MATCH" if d is None else "DIFF", rule, d or "")
    if rule == MZML_INDEXED:
        d = _compare_mzml_stream(pe, pa, exp, act)
        return Result(rel, "MATCH" if d is None else "DIFF", rule, d or "")
    if rule in PHILO_RULES:
        d = PHILO_RULES[rule](pe, pa, exp, act, PurePosixPath(rel).name)
        return Result(rel, "MATCH" if d is None else "DIFF", rule, d or "")
    if rule == FPOP_SET_COLUMN:
        d = compare_fpop_table(pe, pa)
        return Result(rel, "MATCH" if d is None else "DIFF", rule, d or "")
    a = normalize(pe.read_bytes(), rule, exp)
    b = normalize(pa.read_bytes(), rule, act)
    if a == b:
        return Result(rel, "MATCH", rule)
    return Result(rel, "DIFF", rule, _first_diff_bytes(a, b))


def compare_fpop_table(pe: Path, pa: Path) -> str | None:
    """FragPipe_FPOP_Analysis.py:349,420 write ``', '.join(set(peptides))`` into the Peptide Sequence column, so
    the order of that list follows Python's per-process string hashing. Measured on the FPOP golden set:
    12 runs with different PYTHONHASHSEED values; the differences from FragPipe's tables were confined to
    this column and only in the order of the ', '-separated items (FPOP Quant output
    order). That column is compared as a set per row, everything else (header, row count, every other cell)
    byte for byte."""
    rows_e = pe.read_bytes().decode("utf-8").split("\n")
    rows_a = pa.read_bytes().decode("utf-8").split("\n")
    if rows_e[0] != rows_a[0]:
        return f"header differs: expected {rows_e[0]!r} actual {rows_a[0]!r}"
    if len(rows_e) != len(rows_a):
        return f"row count differs: expected {len(rows_e)} actual {len(rows_a)}"
    col = rows_e[0].split("\t").index(FPOP_LIST_COLUMN)
    for n, (re_, ra) in enumerate(zip(rows_e[1:], rows_a[1:]), start=2):
        ce, ca = re_.split("\t"), ra.split("\t")
        if len(ce) != len(ca):
            return f"line {n}: {len(ce)} columns expected, {len(ca)} actual"
        for i, (x, y) in enumerate(zip(ce, ca)):
            same = sorted(x.split(", ")) == sorted(y.split(", ")) if i == col else x == y
            if not same:
                return f"line {n}, column {rows_e[0].split(chr(9))[i]!r}: expected {x[:100]!r} actual {y[:100]!r}"
    return None


def _tree(root: Path) -> dict[str, Path]:
    if not root.is_dir():
        raise FileNotFoundError(f"not a directory: {root}")
    return {p.relative_to(root).as_posix(): p for p in sorted(root.rglob("*")) if p.is_file()}


def check_nonempty(results: list[Result]) -> list[Result]:
    """Comparing nothing is a failure, not a pass."""
    if not any(r.status != "SKIP" for r in results):
        results.append(Result(".", "EMPTY", "-", "no file was compared"))
    return results


def compare_dirs(expected: Path, actual: Path, exp: Side, act: Side, profile: str = "dda") -> list[Result]:
    te, ta = _tree(expected), _tree(actual)
    out = []
    for rel in sorted(set(te) | set(ta)):
        rule = rule_for(rel, profile)
        if rel not in ta:
            out.append(Result(rel, "SKIP" if rule == SKIP else "MISSING", rule, "only in expected"))
        elif rel not in te:
            out.append(Result(rel, "SKIP" if rule == SKIP else "EXTRA", rule, "only in actual"))
        else:
            out.append(compare_files(te[rel], ta[rel], rel, exp, act, profile))
    return check_nonempty(out)


def failed(results: list[Result]) -> bool:
    return any(r.status in FAIL_STATUSES for r in results)


def report(results: list[Result], json_out: Path | None) -> None:
    for r in results:
        print(f"{r.status}\t{r.path}\t{r.rule}\t{r.detail}")
    counts: dict[str, int] = {}
    for r in results:
        counts[r.status] = counts.get(r.status, 0) + 1
    print("SUMMARY\t" + " ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    if json_out:
        json_out.write_text(json.dumps({"results": [asdict(r) for r in results], "counts": counts}, indent=2) + "\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="mode", required=True)
    d = sub.add_parser("dirs", help="compare two directory trees")
    d.add_argument("expected", type=Path)
    d.add_argument("actual", type=Path)
    p = sub.add_parser("pairs", help="compare explicit file pairs")
    p.add_argument("--pair", nargs=3, action="append", metavar=("REL", "EXPECTED", "ACTUAL"), required=True,
                   help="REL is the file's path relative to the run work dir (selects the rule)")
    for s in (d, p):
        s.add_argument("--profile", choices=PROFILES, default="dda", help="rule profile (default dda)")
        s.add_argument("--expected-prefix", action="append", default=[], help="run-specific dir on the expected side (repeatable)")
        s.add_argument("--actual-prefix", action="append", default=[], help="run-specific dir on the actual side (repeatable)")
        s.add_argument("--expected-tmp", action="append", default=[], help="java.io.tmpdir root on the expected side")
        s.add_argument("--actual-tmp", action="append", default=[], help="tmp root on the actual side")
        s.add_argument("--json", type=Path, help="also write results as JSON")
    a = ap.parse_args(argv)
    exp = Side(a.expected_prefix, a.expected_tmp)
    act = Side(a.actual_prefix, a.actual_tmp)
    if a.mode == "dirs":
        results = compare_dirs(a.expected, a.actual, exp, act, a.profile)
    else:
        results = []
        for rel, e, x in a.pair:
            pe, pa = Path(e), Path(x)
            if not pe.is_file():
                raise FileNotFoundError(f"expected file does not exist: {pe}")
            if pa.is_file():
                results.append(compare_files(pe, pa, rel, exp, act, a.profile))
            else:
                results.append(Result(rel, "MISSING", rule_for(rel, a.profile), f"actual file not found: {pa}"))
        results = check_nonempty(results)
    report(results, a.json)
    return 1 if failed(results) else 0


if __name__ == "__main__":
    sys.exit(main())
