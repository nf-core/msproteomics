#!/usr/bin/env python3
"""Comparison of Philosopher filter/report outputs within the measured nondeterminism.

Philosopher 5.1.3-RC9 builds its FDR lists by iterating Go maps, sorts them on
``Probability`` only and takes each tie group's FDR from its first member
(``lib/fil/fdr.go:31-45,132,150``), so q-values, the representative PSM of a
tie and threshold membership change from run to run. The bounds below are the
ones measured over 58 pairs of identical FragPipe 24.0 runs on the same input.
Nothing else is tolerated: every other column and line must match exactly.

Each ``compare_*`` function takes the expected and actual file paths plus the
``Side`` path rewriters and returns ``None`` (match) or a difference description.
Sibling files (``psm.tsv``, ``filter.log``) in the same directory provide the
context the rules need (boundary PSMs, thresholds).
"""

from __future__ import annotations

import bisect
import csv
import re
from pathlib import Path

DOC = "TEST_fragpipe-determinism.md, 'Philosopher filter/report nondeterminism'"
# Row "psm.tsv": |dq| <= 1.01 / max(T_above, 1)
PSM_Q_FACTOR = 1.01
# Row "psm.tsv rows": at most 1 boundary row in either direction
PSM_BOUNDARY_ROWS = 1
# Row "peptide.tsv, ion.tsv Qvalue": one side 0 and other <= FDR cutoff, or |dq| <= 2 / max(T_above, 1)
PEP_ION_Q_FACTOR = 2.0
FDR_CUTOFF = 0.01
# Row "filter.log": counts +-2, decoy +-1, threshold +-2e-4
LOG_COUNT_TOL = 2
LOG_DECOY_TOL = 1
LOG_THRESHOLD_TOL = 2e-4
# Peptide-level converged threshold: the two measured outcomes of PXD001468 commonPTMs_DMO, 0.772 and 0.7725
# (difference 5e-4). Ten reruns of the image's Philosopher v5.1.3-RC9 filter on the golden's inputs gave
# 0.772 (11720 Peptides, decoy=118) 8 times and 0.7725 (11718 Peptides, decoy=117) 2 times; 2000 simulated tie
# orders of the 5.1.3 algorithm gave the same two outcomes (1326 / 674).
PEP_LOG_THRESHOLD_TOL = 5e-4
# Row "peptide.tsv rows": at most 1 boundary peptide in either direction (measured: GGSRGEEVGELSR,
# Probability 0.7720, between the two peptide thresholds above)
PEP_BOUNDARY_ROWS = 1

SPECTRAL_COUNT_COLUMNS = ("Spectral Count", "Total Spectral Count", "Unique Spectral Count", "Razor Spectral Count")
# IonQuant feature columns of peptide/ion/protein.tsv (intensities, retention time and ion mobility of the traced feature)
IONQUANT_FEATURE_RE = re.compile(r"Intensity|Retention|Ion Mobility")
MASS_RE = re.compile(r"\[[^\]]*\]")


class Table:
    """A TSV file: header and rows (lists of str), split on tabs.

    ``ragged``: rows may have another field count than the header, as Philosopher writes its run-level
    tables: every combined_protein.tsv data row ends with a tab the header does not have
    (lib/aba/pro.go: the last field is written as "%s\t" before "\n"), and reprint.*.tsv has the short
    "na" row. Rows are then compared as written, field by field against the header where it has one.
    """

    def __init__(self, path: Path, text: str | None = None, ragged: bool = False):
        text = path.read_text(encoding="utf-8") if text is None else text
        lines = text.split("\n")
        if lines and lines[-1] == "":
            lines = lines[:-1]
        if not lines:
            raise ValueError(f"{path}: empty table")
        self.path = path
        self.header = lines[0].split("\t")
        self.rows = [ln.split("\t") for ln in lines[1:]]
        for n, r in enumerate(self.rows, start=2):
            if not ragged and len(r) < len(self.header):
                # FragPipe's own tables (the psm.tsv PTM-Shepherd rewrites, open search) end rows early: a missing
                # trailing field reads as empty, as in FragPipe's readers
                r.extend([""] * (len(self.header) - len(r)))
            if not ragged and len(r) != len(self.header):
                raise ValueError(f"{path}: line {n} has {len(r)} fields, header has {len(self.header)}")
        self.col = {c: i for i, c in enumerate(self.header)}

    def need(self, *cols: str) -> None:
        missing = [c for c in cols if c not in self.col]
        if missing:
            raise ValueError(f"{self.path}: missing column(s) {missing}")

    def get(self, row: list[str], c: str) -> str:
        return row[self.col[c]]


def _t_above(ref: Table):
    """T_above(row): rows of the reference with Probability strictly greater (doc, rule table preamble).

    Rows without a Probability (IonQuant MBR transfers in a rewritten peptide.tsv) are not ranked.
    """
    probs = sorted(float(v) for v in (ref.get(r, "Probability") for r in ref.rows) if v != "")

    def f(p: float) -> int:
        return len(probs) - bisect.bisect_right(probs, p)

    return f


# ---------------------------------------------------------------- filter.log

COUNT_LINE_RE = re.compile(
    r'^level=info msg="(?:Converged to (?P<fdr>[\d.]+) % FDR with )?(?P<n>\d+) (?P<kind>PSMs|Peptides|Ions)" '
    r"decoy=(?P<decoy>\d+) threshold=(?P<thr>[\d.eE+-]+) total=(?P<total>\d+)$"
)
SUMMARY_LINE_RE = re.compile(r'^level=info msg="(?P<msg>Applying sequential FDR estimation|Final report numbers[^"]*)"(?P<kv>( \w+=\d+)+)$')
PARSING_RE = re.compile(r'^level=info msg="Parsing ')
TIME_RE = re.compile(r'^time="\d{2}:\d{2}:\d{2}" ')


def _log_lines(path: Path, side) -> list[str]:
    text = side.rewrite_paths(path.read_bytes()).decode("utf-8")
    return [TIME_RE.sub("", ln) for ln in text.split("\n")]


def final_psm_threshold(lines: list[str]) -> float | None:
    """Final PSM threshold (doc row "psm.tsv rows"): sequential ``NNNN PSMs`` line, else
    the second ``Converged ... PSMs`` line (2D), else the only one."""
    seq, conv = [], []
    for ln in lines:
        m = COUNT_LINE_RE.match(ln)
        if m and m.group("kind") == "PSMs":
            (conv if m.group("fdr") is not None else seq).append(float(m.group("thr")))
    if seq:
        return seq[-1]
    if len(conv) >= 2:
        return conv[1]
    return conv[0] if conv else None


def final_peptide_threshold(lines: list[str]) -> float | None:
    """Peptide-level converged threshold (the last ``Converged ... Peptides`` line)."""
    thr = [float(m.group("thr")) for m in map(COUNT_LINE_RE.match, lines) if m and m.group("kind") == "Peptides" and m.group("fdr") is not None]
    return thr[-1] if thr else None


def compare_filter_log(pe: Path, pa: Path, exp, act) -> str | None:
    """Doc row "filter.log": paths rewritten, time= stripped, Parsing lines as a sorted set,
    PSM/Peptide/Ion count lines within +-2 / decoy +-1 / threshold +-2e-4, the sequential
    and final summary counts within +-2; every other line (Proteins, charge profiles) exact."""
    le, la = _log_lines(pe, exp), _log_lines(pa, act)
    pars_e = sorted(x for x in le if PARSING_RE.match(x))
    pars_a = sorted(x for x in la if PARSING_RE.match(x))
    if pars_e != pars_a:
        return f"filter.log: Parsing lines differ as a set ({len(pars_e)} vs {len(pars_a)})"
    re_, ra = [x for x in le if not PARSING_RE.match(x)], [x for x in la if not PARSING_RE.match(x)]
    if len(re_) != len(ra):
        return f"filter.log: {len(re_)} vs {len(ra)} non-Parsing lines"
    for n, (x, y) in enumerate(zip(re_, ra), start=1):
        if x == y:
            continue
        mx, my = COUNT_LINE_RE.match(x), COUNT_LINE_RE.match(y)
        if mx and my and mx.group("kind") == my.group("kind") and (mx.group("fdr") is None) == (my.group("fdr") is None):
            nx, ny = int(mx.group("n")), int(my.group("n"))
            dx, dy = int(mx.group("decoy")), int(my.group("decoy"))
            ok = (
                abs(nx - ny) <= LOG_COUNT_TOL
                and abs(dx - dy) <= LOG_DECOY_TOL
                and abs(float(mx.group("thr")) - float(my.group("thr")))
                <= (PEP_LOG_THRESHOLD_TOL if mx.group("kind") == "Peptides" else LOG_THRESHOLD_TOL) + 1e-12
                # total and the FDR percentage are functions of count and decoy
                and int(mx.group("total")) == nx + dx
                and int(my.group("total")) == ny + dy
            )
            if ok:
                continue
            thr_tol = PEP_LOG_THRESHOLD_TOL if mx.group("kind") == "Peptides" else LOG_THRESHOLD_TOL
            return f"filter.log: count line {n} outside +-{LOG_COUNT_TOL}/+-{LOG_DECOY_TOL}/+-{thr_tol}: {x!r} vs {y!r}"
        sx, sy = SUMMARY_LINE_RE.match(x), SUMMARY_LINE_RE.match(y)
        if sx and sy and sx.group("msg") == sy.group("msg"):
            kx = dict(kv.split("=") for kv in sx.group("kv").split())
            ky = dict(kv.split("=") for kv in sy.group("kv").split())
            if kx.keys() == ky.keys() and all(abs(int(kx[k]) - int(ky[k])) <= LOG_COUNT_TOL for k in kx):
                continue
            return f"filter.log: summary line {n} counts outside +-{LOG_COUNT_TOL}: {x!r} vs {y!r}"
        return f"filter.log: line {n} differs: {x!r} vs {y!r}"
    return None


# ------------------------------------------------------------------- psm.tsv


def _psm_keys(t: Table) -> list[tuple[str, ...]]:
    """Doc row "psm.tsv": key on Spectrum, plus Spectrum File when Spectrum is not unique."""
    t.need("Spectrum")
    spectra = [t.get(r, "Spectrum") for r in t.rows]
    if len(set(spectra)) == len(spectra):
        return [(s,) for s in spectra]
    t.need("Spectrum File")
    return [(t.get(r, "Spectrum"), t.get(r, "Spectrum File")) for r in t.rows]


def _sibling_threshold(path: Path, side) -> float | None:
    log = path.parent / "filter.log"
    return final_psm_threshold(_log_lines(log, side)) if log.is_file() else None


def boundary_psms(pe: Path, pa: Path, exp, act) -> tuple[list[list[str]], list[list[str]], Table, Table]:
    """Rows only in expected / only in actual psm.tsv (the boundary PSMs)."""
    te, ta = Table(pe), Table(pa)
    ke, ka = _psm_keys(te), _psm_keys(ta)
    se, sa = set(ke), set(ka)
    only_e = [r for k, r in zip(ke, te.rows) if k not in sa]
    only_a = [r for k, r in zip(ka, ta.rows) if k not in se]
    return only_e, only_a, te, ta


def compare_psm(pe: Path, pa: Path, exp, act) -> str | None:
    """Doc rows "psm.tsv" and "psm.tsv rows"."""
    only_e, only_a, te, ta = boundary_psms(pe, pa, exp, act)
    if te.header != ta.header:
        return "psm.tsv: header differs"
    te.need("Probability", "Qvalue")
    if only_e or only_a:
        if len(only_e) > PSM_BOUNDARY_ROWS or len(only_a) > PSM_BOUNDARY_ROWS:
            return f"psm.tsv: {len(only_e)} rows only in expected, {len(only_a)} only in actual (max {PSM_BOUNDARY_ROWS} each)"
        t1, t2 = _sibling_threshold(pe, exp), _sibling_threshold(pa, act)
        if t1 is None or t2 is None:
            return "psm.tsv: rows differ and no sibling filter.log threshold to check them against"
        lo, hi = min(t1, t2), max(t1, t2)
        for r, t in [(r, te) for r in only_e] + [(r, ta) for r in only_a]:
            p = float(t.get(r, "Probability"))
            if not lo <= p <= hi:
                return f"psm.tsv: boundary row {t.get(r, 'Spectrum')} Probability {p} not within thresholds {lo}..{hi}"
    ke, ka = _psm_keys(te), _psm_keys(ta)
    common = set(ke) & set(ka)
    re_ = [(k, r) for k, r in zip(ke, te.rows) if k in common]
    ra = [(k, r) for k, r in zip(ka, ta.rows) if k in common]
    t_above = _t_above(te)
    qi = te.col["Qvalue"]
    for (k1, r1), (k2, r2) in zip(re_, ra):
        if k1 != k2:
            return f"psm.tsv: row order differs at {k1} vs {k2}"
        for i, (x, y) in enumerate(zip(r1, r2)):
            if x == y:
                continue
            if i != qi:
                return f"psm.tsv: {k1[0]} column {te.header[i]!r}: {x!r} vs {y!r}"
            bound = PSM_Q_FACTOR / max(t_above(float(te.get(r1, "Probability"))), 1)
            if abs(float(x) - float(y)) > bound:
                return f"psm.tsv: {k1[0]} Qvalue {x} vs {y} exceeds {PSM_Q_FACTOR}/max(T_above,1) = {bound:.3g}"
    return None


# --------------------------------------------- peptide / ion / protein tables


def _stripped(t: Table, r: list[str]) -> str:
    """Doc row "ion.tsv": Modified Sequence without [mass] and n[...], or Peptide Sequence."""
    mod = t.get(r, "Modified Sequence")
    return MASS_RE.sub("", mod).lstrip("n") if mod else t.get(r, "Peptide Sequence")


def _mod_signature(t: Table, r: list[str]) -> str:
    mod = t.get(r, "Modified Sequence")
    if not mod:
        return ""
    return _stripped(t, r) + "|" + ",".join(sorted(MASS_RE.findall(mod)))


def _boundary_counts(pe: Path, pa: Path, exp, act, key) -> dict:
    """Doc row "Spectral counts": boundary PSMs (sibling psm.tsv symmetric difference) per entry key."""
    se, sa = pe.parent / "psm.tsv", pa.parent / "psm.tsv"
    if not (se.is_file() and sa.is_file()):
        return {}
    only_e, only_a, te, ta = boundary_psms(se, sa, exp, act)
    counts: dict = {}
    for r, t in [(r, te) for r in only_e] + [(r, ta) for r in only_a]:
        for k in key(t, r):
            counts[k] = counts.get(k, 0) + 1
    return counts


def _psm_peptide_keys(t: Table, r: list[str]) -> list:
    return [t.get(r, "Peptide")]


def _psm_ion_keys(t: Table, r: list[str]) -> list:
    return [(t.get(r, "Peptide"), t.get(r, "Charge"))]


def _psm_protein_keys(t: Table, r: list[str]) -> list:
    keys = [t.get(r, "Protein")]
    if "Mapped Proteins" in t.col:
        keys += [x.strip() for x in t.get(r, "Mapped Proteins").split(",") if x.strip()]
    return keys


def _psm_ion_evidence(table_path: Path) -> dict:
    """Per ion (Peptide, Charge) of the sibling psm.tsv: its PSMs' Modified Peptide and Expectation values."""
    sp = table_path.parent / "psm.tsv"
    if not sp.is_file():
        return {}
    t = Table(sp)
    t.need("Peptide", "Charge", "Modified Peptide", "Expectation")
    ev: dict = {}
    for r in t.rows:
        e = ev.setdefault((t.get(r, "Peptide"), t.get(r, "Charge")), {"Modified Sequence": set(), "Expectation": set()})
        e["Modified Sequence"].add(t.get(r, "Modified Peptide") or t.get(r, "Peptide"))
        e["Expectation"].add(float(t.get(r, "Expectation")))
    return ev


def _representative_of_ion(column: str, value: str, t: Table, row: list[str], evidence: dict) -> bool:
    """Whether an ion.tsv Modified Sequence / Expectation is the value of one of the ion's own PSMs (Philosopher
    picks the representative PSM run to run). An ion without PSMs (Spectral Count 0) has an Expectation from
    PSMs that did not pass the filter: measured to vary too, with nothing to check it against."""
    if column == "Expectation" and t.get(row, "Spectral Count") == "0":
        return True
    e = evidence.get((t.get(row, "Peptide Sequence"), t.get(row, "Charge")))
    if e is None:
        return False
    return (float(value) if column == "Expectation" else value) in e[column]


def _drop_boundary_peptides(pe: Path, pa: Path, exp, act, te: Table, ta: Table):
    """Doc row "peptide.tsv rows": the peptide-level FDR threshold moves between runs, so a peptide whose
    Probability lies between the two sibling filter.log thresholds may be present on one side only
    (at most PEP_BOUNDARY_ROWS per side). Returns (expected rows, actual rows) without them, or an error."""
    te.need("Peptide", "Probability")
    ke, ka = [r[te.col["Peptide"]] for r in te.rows], [r[ta.col["Peptide"]] for r in ta.rows]
    se, sa = set(ke), set(ka)
    only_e = [r for k, r in zip(ke, te.rows) if k not in sa]
    only_a = [r for k, r in zip(ka, ta.rows) if k not in se]
    if len(only_e) > PEP_BOUNDARY_ROWS or len(only_a) > PEP_BOUNDARY_ROWS:
        return None, None, f"peptide.tsv: {len(only_e)} rows only in expected, {len(only_a)} only in actual (max {PEP_BOUNDARY_ROWS} each)"
    logs = [pe.parent / "filter.log", pa.parent / "filter.log"]
    if not all(f.is_file() for f in logs):
        return None, None, "peptide.tsv: rows differ and no sibling filter.log threshold to check them against"
    t1, t2 = (final_peptide_threshold(_log_lines(f, side)) for f, side in zip(logs, (exp, act)))
    if t1 is None or t2 is None:
        return None, None, "peptide.tsv: rows differ and filter.log has no peptide-level converged threshold"
    lo, hi = min(t1, t2), max(t1, t2)
    # The side with the lower threshold keeps more peptides: only that side may hold the extra row.
    if (only_e or only_a) and t1 == t2:
        return None, None, f"peptide.tsv: rows only on one side but both peptide-level thresholds are {t1}"
    extra_side = only_e if t1 < t2 else only_a
    other_side = only_a if t1 < t2 else only_e
    if other_side:
        return None, None, f"peptide.tsv: row(s) only on the side with the higher peptide-level threshold ({hi}; lower {lo})"
    for r in extra_side:
        t = te if t1 < t2 else ta
        p = float(t.get(r, "Probability"))
        if not lo <= p <= hi:
            return None, None, f"peptide.tsv: boundary row {t.get(r, 'Peptide')} Probability {p} not within thresholds {lo}..{hi}"
    rows_e = [r for k, r in zip(ke, te.rows) if k in sa]
    rows_a = [r for k, r in zip(ka, ta.rows) if k in se]
    if len(rows_e) != len(rows_a):
        return None, None, f"peptide.tsv: {len(rows_e)} vs {len(rows_a)} rows after dropping boundary peptides"
    return rows_e, rows_a, None


def compare_entry_table(pe: Path, pa: Path, exp, act, kind: str, full: bool) -> str | None:
    """peptide.tsv / ion.tsv / protein.tsv.

    Always (doc row "Spectral counts"; dda and dia): spectral-count columns may differ
    by the number of boundary PSMs mapped to the entry; ``protein.tsv`` otherwise exact.
    peptide.tsv, every scope: Prev AA/Next AA dropped and the Qvalue rule (measured
    Philosopher run-to-run variation, see the peptide.tsv branch).
    ion.tsv, every scope: keyed on (stripped sequence, Charge, M/Z), rows sorted, Prev AA/Next AA
    dropped, the Qvalue rule, and Modified Sequence / Expectation each a value of one of the ion's
    own PSMs (sibling psm.tsv) on its side (measured, see the ion.tsv branch).
    ``full`` (dia only, doc "Workflow scope"): ion.tsv Expectation dropped and Modified Sequence also
    accepted as an equal mass multiset.
    """
    te, ta = Table(pe), Table(pa)
    if te.header != ta.header:
        return f"{kind}: header differs"
    drop: set[str] = set()
    qcol = None
    if kind == "peptide.tsv":
        te.need("Peptide")
        keyf = lambda t, r: (t.get(r, "Peptide"),)  # noqa: E731
        bkey, bmap = _psm_peptide_keys, (lambda t, r: t.get(r, "Peptide"))
        # Prev AA/Next AA and Qvalue in every scope: Philosopher filter+report rerun on identical
        # inputs (PXD001468 Basic-Search, 6 runs) gave LIEALLR Next AA "A" in 5 runs and "L" in 1,
        # and on commonPTMs_DMO (8 runs) AANAIQQQTVVVPASSLANAK Qvalue 0.00123903926801 and
        # 0.00123974823574 in 4 runs each. IonQuant rewrites the intensity columns, not these.
        drop = {"Prev AA", "Next AA"}
        qcol = "Qvalue"
    elif kind == "ion.tsv":
        te.need("Charge", "M/Z", "Modified Sequence", "Peptide Sequence")
        keyf = lambda t, r: (_stripped(t, r), t.get(r, "Charge"), t.get(r, "M/Z"))  # noqa: E731
        bkey, bmap = _psm_ion_keys, (lambda t, r: (_stripped(t, r), t.get(r, "Charge")))
        # Every scope: Philosopher report on identical inputs (PXD001468 commonPTMs_DMO 6 runs, Basic-Search
        # 4 runs, FragPipe 24.0 image) writes ion.tsv in a different row order each run, with another
        # representative PSM per ion: Expectation (~1400 rows), Qvalue (~1000), Modified Sequence (~100
        # with --mapmods, ~25 without) and Prev/Next AA (1-3) differ; the ion keys never do.
        drop = {"Expectation", "Prev AA", "Next AA"} if full else {"Prev AA", "Next AA"}
        qcol = "Qvalue"
    elif kind == "protein.tsv":
        te.need("Protein")
        keyf = lambda t, r: (t.get(r, "Protein"),)  # noqa: E731
        bkey, bmap = _psm_protein_keys, None
    else:
        raise ValueError(kind)
    keyed = False
    if kind == "ion.tsv":
        # Doc row "ion.tsv": key on (stripped sequence, Charge, M/Z) and sort rows. Philosopher's ion.tsv has
        # one row per key; an IonQuant-rewritten ion.tsv lists localisations of one key separately, in a
        # fixed order, and is compared row by row.
        ke = [keyf(te, r) for r in te.rows]
        ka = [keyf(ta, r) for r in ta.rows]
        keyed = len(set(ke)) == len(ke) and len(set(ka)) == len(ka)
        if full and not keyed:
            raise ValueError(f"{pe}: ion key (stripped sequence, Charge, M/Z) is not unique")
    if keyed:
        if sorted(ke) != sorted(ka):
            return f"{kind}: row keys differ ({len(ke)} vs {len(ka)} rows)"
        rows_a = dict(zip(ka, ta.rows))
        pairs = [(k, r, rows_a[k]) for k, r in zip(ke, te.rows)]
    else:
        # Rows in the same order; the row key is only used in messages.
        if len(te.rows) != len(ta.rows):
            if kind != "peptide.tsv":
                return f"{kind}: {len(te.rows)} vs {len(ta.rows)} rows"
            te.rows, ta.rows, err = _drop_boundary_peptides(pe, pa, exp, act, te, ta)
            if err:
                return err
        pairs = [(keyf(te, r1), r1, r2) for r1, r2 in zip(te.rows, ta.rows)]
    sc_cols = {c for c in SPECTRAL_COUNT_COLUMNS if c in te.col}
    boundary = _boundary_counts(pe, pa, exp, act, bkey) if sc_cols else {}
    t_above = _t_above(te) if qcol else None
    ions_e, ions_a = (_psm_ion_evidence(pe), _psm_ion_evidence(pa)) if kind == "ion.tsv" else ({}, {})
    for k, r1, r2 in pairs:
        for i, (x, y) in enumerate(zip(r1, r2)):
            c = te.header[i]
            if x == y or c in drop:
                continue
            if kind == "ion.tsv" and c in ("Modified Sequence", "Expectation"):
                if c == "Modified Sequence" and full and _mod_signature(te, r1) == _mod_signature(ta, r2):
                    continue
                if _representative_of_ion(c, x, te, r1, ions_e) and _representative_of_ion(c, y, ta, r2, ions_a):
                    continue
                return f"{kind}: {k} {c} {x!r} vs {y!r} is not a value of the ion's own PSMs"
            allowed = boundary.get(k[0], 0) if kind == "protein.tsv" else boundary.get(bmap(te, r1), 0)
            if allowed and IONQUANT_FEATURE_RE.search(c):
                # IonQuant quantifies the entry from its PSMs: a boundary PSM changes its feature values
                # (measured on a 12-run LFQ set, experiment_1_1: FragPipe's filter on identical inputs gives 15014 or
                # 15016 PSMs, 5 vs 7 of 12 runs; the extra PSM changes only its peptide's, ion's and
                # protein's Intensity / Retention Time FWHM / Ion Mobility End)
                continue
            if c in sc_cols:
                if abs(int(x) - int(y)) <= allowed:
                    continue
                return f"{kind}: {k} {c} {x} vs {y}; {allowed} boundary PSM(s) map to it"
            if c == qcol:
                if "" in (x, y, te.get(r1, "Probability")):
                    # no Philosopher probability (MBR row): no envelope, the values must be equal
                    return f"{kind}: {k} Qvalue {x!r} vs {y!r} (row without a Philosopher Probability)"
                qx, qy = float(x), float(y)
                if (qx == 0 and qy <= FDR_CUTOFF) or (qy == 0 and qx <= FDR_CUTOFF):
                    continue
                bound = PEP_ION_Q_FACTOR / max(t_above(float(te.get(r1, "Probability"))), 1)
                if abs(qx - qy) <= bound:
                    continue
                return f"{kind}: {k} Qvalue {x} vs {y} exceeds {PEP_ION_Q_FACTOR}/max(T_above,1) = {bound:.3g}"
            return f"{kind}: {k} column {c!r}: {x!r} vs {y!r}"
        if not keyed and keyf(te, r1) != keyf(ta, r2):
            return f"{kind}: row order differs at {k}"
    return None


# ---------------------------------------------------------------- run-level tables

RUN_TABLES = (
    "combined_protein.tsv", "combined_peptide.tsv", "combined_modified_peptide.tsv", "combined_ion.tsv",
    "msstats.csv", "msstats_ptm.csv", "reprint.spc.tsv", "reprint.int.tsv",
)
# The columns a boundary PSM changes (measured on a 12-run LFQ set, experiment_1_1, FragPipe's filter on identical
# inputs 15014 vs 15016 PSMs in 5 vs 7 of 12 runs): spectral counts, IonQuant intensities, RT FWHM, ion mobility.
RUN_TABLE_VALUE_RE = re.compile(r"Spectral Count|Intensity|Retention Time FWHM|Ion Mobility|_SPC$|_INT$")


def _rows(path: Path) -> tuple[list[str], list[list[str]]]:
    if path.suffix == ".csv":
        with path.open(encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh))
    else:
        t = Table(path, ragged=True)
        rows = [t.header] + t.rows
    if not rows:
        raise ValueError(f"{path}: empty table")
    return rows[0], rows[1:]


def run_boundary_psms(pe: Path, pa: Path, exp, act) -> tuple[list[dict], str | None]:
    """Boundary PSMs of every group dir next to a run-level table (both sides' <group>/psm.tsv).

    Each group's psm.tsv must pass compare_psm (its boundary rows are within the measured envelope).
    """
    found = []
    for se in sorted(pe.parent.glob("*/psm.tsv")):
        group = se.parent.name
        sa = pa.parent / group / "psm.tsv"
        if not sa.is_file():
            return [], f"{group}/psm.tsv missing on the actual side"
        err = compare_psm(se, sa, exp, act)
        if err:
            return [], f"{group}/{err}"
        only_e, only_a, te, ta = boundary_psms(se, sa, exp, act)
        for r, t in [(r, te) for r in only_e] + [(r, ta) for r in only_a]:
            proteins = {t.get(r, "Protein")}
            if "Mapped Proteins" in t.col:
                proteins |= {x.strip() for x in t.get(r, "Mapped Proteins").split(",") if x.strip()}
            found.append({
                "group": group,
                "run": t.get(r, "Spectrum").split(".")[0],
                "peptide": t.get(r, "Peptide"),
                "charge": t.get(r, "Charge"),
                "proteins": proteins,
                "protein_ids": {p.split("|")[1] for p in proteins if p.count("|") >= 2},
            })
    return found, None


def _run_row_entry(name: str, header: list[str], row: list[str], b: dict) -> bool:
    """Whether a run-level table row is an entry boundary PSM ``b`` maps to."""
    get = dict(zip(header, row)).get
    if name == "combined_protein.tsv":
        return get("Protein") in b["proteins"]
    if name in ("combined_peptide.tsv", "combined_modified_peptide.tsv"):
        return get("Peptide Sequence") == b["peptide"]
    if name == "combined_ion.tsv":
        return get("Peptide Sequence") == b["peptide"] and get("Charge") == b["charge"]
    if name in ("msstats.csv", "msstats_ptm.csv"):
        return (MASS_RE.sub("", get("PeptideSequence") or "") == b["peptide"] and get("PrecursorCharge") == b["charge"]
                and get("Run") == b["run"] and get("ProteinName") in b["proteins"])
    if name in ("reprint.spc.tsv", "reprint.int.tsv"):
        return get("PROTID") in b["protein_ids"]
    raise ValueError(name)


def _run_column_of(name: str, column: str, b: dict) -> bool:
    """Whether a differing column is one boundary PSM ``b`` changes (its group's value columns)."""
    if not RUN_TABLE_VALUE_RE.search(column):
        return False
    if name in ("msstats.csv", "msstats_ptm.csv"):
        return column == "Intensity"  # one row per run, matched on Run
    return column.startswith((f"{b['group']} ", f"{b['group']}_", "Combined "))


def compare_run_table(pe: Path, pa: Path, exp, act, name: str) -> str | None:
    """combined_*.tsv, msstats*.csv, reprint.*.tsv: byte-exact rows, except that a row of an entry a boundary
    PSM maps to may differ in that PSM's group's spectral-count / intensity / RT FWHM / ion-mobility columns."""
    he, re_ = _rows(pe)
    ha, ra = _rows(pa)
    if he != ha:
        return f"{name}: header differs"
    if len(re_) != len(ra):
        return f"{name}: {len(re_)} vs {len(ra)} rows"
    boundary = None
    for n, (r1, r2) in enumerate(zip(re_, ra), start=2):
        if r1 == r2:
            continue
        if boundary is None:
            boundary, err = run_boundary_psms(pe, pa, exp, act)
            if err:
                return f"{name}: {err}"
        # fields past the header (Philosopher's trailing tab) have no column name and are never explained
        width = max(len(r1), len(r2))
        x1, x2 = r1 + [None] * (width - len(r1)), r2 + [None] * (width - len(r2))
        cols = [he[i] if i < len(he) else f"field {i + 1}" for i in range(width) if x1[i] != x2[i]]
        if not any(_run_row_entry(name, he, r1, b) and all(_run_column_of(name, c, b) for c in cols) for b in boundary):
            return f"{name}: line {n} ({r1[0]}) differs in {cols[:6]} and no boundary PSM explains it"
    return None
