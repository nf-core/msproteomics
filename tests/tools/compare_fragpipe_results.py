#!/usr/bin/env python3
"""Compare an msproteomics result directory against a FragPipe result directory.

Both directories are searched recursively for FragPipe-style outputs:

* per-experiment ``psm.tsv`` / ``peptide.tsv`` / ``ion.tsv`` / ``protein.tsv``
  (the experiment name is the parent directory name),
* ``combined_protein.tsv`` / ``combined_peptide.tsv`` /
  ``combined_modified_peptide.tsv`` / ``combined_ion.tsv`` (IonQuant / FragPipe),
* TMT-Integrator tables ``abundance_*.tsv`` / ``ratio_*.tsv`` (paired by file name),
* ``fragger.params`` / ``fragpipe.workflow`` (diffed when present on both sides).

Writes ``comparison.md`` and ``comparison.json`` to ``--outdir`` and exits with
0 (PASS), 1 (FAIL) or 2 (input error).

Stdlib only so it runs in any container with Python >= 3.8.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
import time
from pathlib import Path

csv.field_size_limit(sys.maxsize)

EXIT_PASS = 0
EXIT_FAIL = 1
EXIT_ERROR = 2

EXPERIMENT_FILES = ("psm.tsv", "peptide.tsv", "ion.tsv", "protein.tsv")
COMBINED_FILES = (
    "combined_protein.tsv",
    "combined_peptide.tsv",
    "combined_modified_peptide.tsv",
    "combined_ion.tsv",
)
TMT_PATTERN = re.compile(r"^(abundance|ratio)_.*\.tsv$")
PARAM_FILES = ("fragger.params", "fragpipe.workflow")
# Combined-table per-sample quant columns are "<experiment> <suffix>".
# Longest suffix first so "MaxLFQ Intensity" is not taken as "Intensity".
COMBINED_QUANT_SUFFIXES = (" MaxLFQ Intensity", " Intensity")
# TMT-Integrator writes channel columns immediately after this column.
TMT_LAST_META_COLUMN = "ReferenceIntensity"
# Parameter keys that legitimately differ between environments (paths, threads).
ENV_PARAM_PATTERN = re.compile(
    r"(database_name|db-path|num_threads|^threads$|\.ram$|^ram$|workdir|output[-_]?(dir|location)|\.dir$|-dir$)",
    re.IGNORECASE,
)
ASSIGNED_MOD_PATTERN = re.compile(r"^(\d*)([A-Za-z]|N-term|C-term)\((-?[0-9.]+)\)$")


class InputError(Exception):
    """Raised for missing files, missing columns or unmappable inputs."""


# --------------------------------------------------------------------------- #
# Reading
# --------------------------------------------------------------------------- #


def read_tsv(path: Path, required: list[str]) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV and verify that the required columns exist."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        if header is None:
            raise InputError(f"{path}: file is empty (no header)")
        missing = [c for c in required if c not in header]
        if missing:
            raise InputError(f"{path}: missing required column(s) {missing}")
        return list(header), list(reader)


def not_decoy(row: dict[str, str]) -> bool:
    return row.get("Is Decoy", "false").strip().lower() != "true"


def mod_key(sequence: str, assigned: str, decimals: int, path: Path) -> str:
    """Canonical modified-peptide key from 'Assigned Modifications'.

    'Modified Peptide' notation differs between tables ([147] residue+mod mass in
    psm.tsv, [57.0215] delta mass in ion/combined tables), so the key is built
    from 'Assigned Modifications' (e.g. '3M(15.9949), N-term(42.0106)').
    """
    mods = []
    for token in (t.strip() for t in assigned.split(",")):
        if not token:
            continue
        match = ASSIGNED_MOD_PATTERN.match(token)
        if match is None:
            raise InputError(f"{path}: cannot parse assigned modification {token!r}")
        pos, residue, mass = match.groups()
        if residue == "N-term":
            pos = "0"
        elif residue == "C-term":
            pos = str(len(sequence) + 1)
        mods.append((int(pos), residue, f"{float(mass):.{decimals}f}"))
    mods.sort()
    return sequence + "".join(f"|{p}{r}({m})" for p, r, m in mods)


def to_float(value: str) -> float | None:
    value = value.strip()
    if not value:
        return None
    return float(value)


def discover(root: Path, exclude: re.Pattern | None) -> dict:
    """Locate FragPipe-style outputs below ``root``."""
    if not root.is_dir():
        raise InputError(f"{root}: not a directory")
    experiments: dict[str, dict[str, Path]] = {}
    combined: dict[str, Path] = {}
    tmt: dict[str, Path] = {}
    params: dict[str, Path] = {}

    def claim(store: dict, key: str, path: Path, what: str) -> None:
        if key in store:
            raise InputError(
                f"{root}: duplicate {what} {key!r}: {store[key]} and {path} "
                "(use --exclude to drop one)"
            )
        store[key] = path

    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        dirnames.sort()
        for name in sorted(filenames):
            path = Path(dirpath) / name
            if exclude is not None and exclude.search(str(path)):
                continue
            if name in EXPERIMENT_FILES:
                exp = experiments.setdefault(Path(dirpath).name, {})
                claim(exp, name, path, f"experiment file in {Path(dirpath).name}")
            elif name in COMBINED_FILES:
                claim(combined, name, path, "combined table")
            elif TMT_PATTERN.match(name):
                claim(tmt, name, path, "TMT-Integrator table")
            elif name in PARAM_FILES:
                claim(params, name, path, "parameter file")
    experiments = {k: v for k, v in experiments.items() if "psm.tsv" in v}
    if not experiments and not combined and not tmt:
        raise InputError(f"{root}: no psm.tsv, combined_*.tsv or TMT-Integrator tables found")
    return {"experiments": experiments, "combined": combined, "tmt": tmt, "params": params}


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #


def pearson(x: list[float], y: list[float]) -> float | None:
    n = len(x)
    if n < 2:
        return None
    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        # Identical constant vectors agree perfectly; otherwise undefined.
        return 1.0 if x == y else None
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sxx * syy)


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            result[order[k]] = avg
        i = j + 1
    return result


def spearman(x: list[float], y: list[float]) -> float | None:
    if len(x) < 2:
        return None
    return pearson(ranks(x), ranks(y))


def median(values: list[float]) -> float | None:
    if not values:
        return None
    s = sorted(values)
    mid = len(s) // 2
    return s[mid] if len(s) % 2 else (s[mid - 1] + s[mid]) / 2


def set_overlap(a: set, b: set) -> dict:
    union = a | b
    return {
        "test": len(a),
        "reference": len(b),
        "shared": len(a & b),
        "only_test": len(a - b),
        "only_reference": len(b - a),
        "jaccard": (len(a & b) / len(union)) if union else 1.0,
    }


def quant_agreement(pairs: list[tuple[float | None, float | None]], min_n: int, log: bool) -> dict:
    """Agreement of paired intensities; zero/blank values are treated as missing."""
    xs: list[float] = []
    ys: list[float] = []
    only_test = only_ref = 0
    for a, b in pairs:
        a_ok = a is not None and a > 0 if log else a is not None
        b_ok = b is not None and b > 0 if log else b is not None
        if a_ok and b_ok:
            xs.append(math.log2(a) if log else a)
            ys.append(math.log2(b) if log else b)
        elif a_ok:
            only_test += 1
        elif b_ok:
            only_ref += 1
    diffs = [a - b for a, b in zip(xs, ys)]
    result = {
        "n": len(xs),
        "only_test": only_test,
        "only_reference": only_ref,
        "pearson": None,
        "spearman": None,
        "median_log2_ratio": median(diffs),
        "pct_within_2fold": (100.0 * sum(abs(d) <= 1 for d in diffs) / len(diffs)) if diffs else None,
        "note": None,
    }
    if len(xs) < min_n:
        result["note"] = f"n={len(xs)} < --min-quant-n={min_n}; correlation not computed"
    else:
        result["pearson"] = pearson(xs, ys)
        result["spearman"] = spearman(xs, ys)
    return result


# --------------------------------------------------------------------------- #
# Table loaders: each returns {id_kind: set} and {quant_name: {key: value}}
# --------------------------------------------------------------------------- #


def load_psm(path: Path, decimals: int) -> dict:
    _, rows = read_tsv(path, ["Spectrum", "Peptide", "Assigned Modifications", "Protein ID"])
    rows = [r for r in rows if not_decoy(r)]
    by_spectrum = {}
    for r in rows:
        if r["Spectrum"] in by_spectrum:
            raise InputError(f"{path}: duplicate Spectrum {r['Spectrum']!r}")
        by_spectrum[r["Spectrum"]] = r
    return {
        "rows": by_spectrum,
        "ids": {
            "peptide_sequences": {r["Peptide"] for r in rows},
            "modified_peptides": {mod_key(r["Peptide"], r["Assigned Modifications"], decimals, path) for r in rows},
            "proteins": {r["Protein ID"] for r in rows},
        },
        "count": len(rows),
    }


def load_keyed(path: Path, key_cols: list[str], value_cols: list[str], decimals: int, mod_cols: tuple[str, str] | None) -> dict:
    """Load a table keyed by key_cols (optionally canonicalising mods)."""
    required = list(key_cols) + (list(mod_cols) if mod_cols else [])
    header, rows = read_tsv(path, required)
    rows = [r for r in rows if not_decoy(r)]
    present = [c for c in value_cols if c in header]
    keys = {}
    for r in rows:
        parts = [r[c] for c in key_cols]
        if mod_cols:
            parts.insert(0, mod_key(r[mod_cols[0]], r[mod_cols[1]], decimals, path))
        key = "/".join(parts)
        if key in keys:
            raise InputError(f"{path}: duplicate key {key!r}")
        keys[key] = {c: to_float(r[c]) for c in present}
    return {"keys": keys, "value_cols": present, "header": header}


EXPERIMENT_TABLES = {
    # file: (key columns, quant columns, mod columns or None, id kind)
    "peptide.tsv": (["Peptide"], ["Intensity"], None, "peptide_sequences"),
    "ion.tsv": (["Charge"], ["Intensity"], ("Peptide Sequence", "Assigned Modifications"), "ions"),
    "protein.tsv": (["Protein ID"], ["Total Intensity", "Unique Intensity", "Razor Intensity"], None, "proteins"),
}
COMBINED_TABLES = {
    "combined_protein.tsv": (["Protein ID"], None, "proteins"),
    "combined_peptide.tsv": (["Peptide Sequence"], None, "peptide_sequences"),
    "combined_modified_peptide.tsv": ([], ("Peptide Sequence", "Assigned Modifications"), "modified_peptides"),
    "combined_ion.tsv": (["Charge"], ("Peptide Sequence", "Assigned Modifications"), "ions"),
}


def compare_values(a: dict, b: dict, cols: list[str], min_n: int, log: bool) -> dict:
    out = {}
    shared = a["keys"].keys() & b["keys"].keys()
    for col in cols:
        pairs = [(a["keys"][k][col], b["keys"][k][col]) for k in shared]
        out[col] = quant_agreement(pairs, min_n, log)
    return out


def quant_cols_note(path_a: Path, path_b: Path, cols_a: list[str], cols_b: list[str], wanted: list[str]) -> list[str]:
    absent = [c for c in wanted if c not in cols_a and c not in cols_b]
    one_sided = [c for c in wanted if (c in cols_a) != (c in cols_b)]
    if one_sided:
        raise InputError(f"quant column(s) {one_sided} present in only one of {path_a} / {path_b}")
    return absent


def compare_experiment(name_t: str, name_r: str, files_t: dict, files_r: dict, args) -> dict:
    result = {"test_name": name_t, "reference_name": name_r, "counts": {}, "overlap": {}, "quant": {}, "absent_quant_columns": {}}
    psm_t = load_psm(files_t["psm.tsv"], args.mod_decimals)
    psm_r = load_psm(files_r["psm.tsv"], args.mod_decimals)
    result["counts"]["psms"] = (psm_t["count"], psm_r["count"])
    for kind in ("peptide_sequences", "modified_peptides", "proteins"):
        result["overlap"][f"psm:{kind}"] = set_overlap(psm_t["ids"][kind], psm_r["ids"][kind])
    result["counts"]["modified_peptides"] = (len(psm_t["ids"]["modified_peptides"]), len(psm_r["ids"]["modified_peptides"]))

    shared = psm_t["rows"].keys() & psm_r["rows"].keys()
    same_pep = sum(psm_t["rows"][s]["Peptide"] == psm_r["rows"][s]["Peptide"] for s in shared)
    same_mod = sum(
        mod_key(psm_t["rows"][s]["Peptide"], psm_t["rows"][s]["Assigned Modifications"], args.mod_decimals, files_t["psm.tsv"])
        == mod_key(psm_r["rows"][s]["Peptide"], psm_r["rows"][s]["Assigned Modifications"], args.mod_decimals, files_r["psm.tsv"])
        for s in shared
    )
    psm_cmp = {
        "spectra": set_overlap(set(psm_t["rows"]), set(psm_r["rows"])),
        "pct_same_peptide": (100.0 * same_pep / len(shared)) if shared else None,
        "pct_same_modified_peptide": (100.0 * same_mod / len(shared)) if shared else None,
        "scores": {},
    }
    for col in ("Probability", "Hyperscore", "Expectation"):
        in_t = col in next(iter(psm_t["rows"].values()), {})
        in_r = col in next(iter(psm_r["rows"].values()), {})
        if in_t and in_r:
            same = [s for s in shared if psm_t["rows"][s]["Peptide"] == psm_r["rows"][s]["Peptide"]]
            pairs = [(to_float(psm_t["rows"][s][col]), to_float(psm_r["rows"][s][col])) for s in same]
            psm_cmp["scores"][col] = quant_agreement(pairs, args.min_quant_n, log=False)
    result["psm"] = psm_cmp
    pairs = [(to_float(psm_t["rows"][s].get("Intensity", "")), to_float(psm_r["rows"][s].get("Intensity", ""))) for s in shared]
    if "Intensity" in next(iter(psm_t["rows"].values()), {}) and "Intensity" in next(iter(psm_r["rows"].values()), {}):
        result["quant"]["psm.tsv:Intensity"] = quant_agreement(pairs, args.min_quant_n, log=True)

    for fname, (key_cols, value_cols, mod_cols, kind) in EXPERIMENT_TABLES.items():
        if (fname in files_t) != (fname in files_r):
            raise InputError(f"{fname} present for only one side of experiment {name_t}/{name_r}: {files_t.get(fname) or files_r.get(fname)}")
        if fname not in files_t:
            continue
        a = load_keyed(files_t[fname], key_cols, value_cols, args.mod_decimals, mod_cols)
        b = load_keyed(files_r[fname], key_cols, value_cols, args.mod_decimals, mod_cols)
        count_name = fname.replace(".tsv", "s")
        result["counts"][count_name] = (len(a["keys"]), len(b["keys"]))
        result["overlap"][f"{fname}:{kind}"] = set_overlap(set(a["keys"]), set(b["keys"]))
        result["absent_quant_columns"][fname] = quant_cols_note(files_t[fname], files_r[fname], a["value_cols"], b["value_cols"], value_cols)
        for col, metrics in compare_values(a, b, a["value_cols"], args.min_quant_n, log=True).items():
            result["quant"][f"{fname}:{col}"] = metrics
    return result


def combined_samples(header: list[str]) -> dict[str, dict[str, str]]:
    """Map quant suffix -> {experiment name: column}."""
    out: dict[str, dict[str, str]] = {}
    for col in header:
        for suffix in COMBINED_QUANT_SUFFIXES:
            if col.endswith(suffix):
                out.setdefault(suffix.strip(), {})[col[: -len(suffix)]] = col
                break
    return out


def compare_combined(fname: str, path_t: Path, path_r: Path, exp_map: dict[str, str], args) -> dict:
    key_cols, mod_cols, kind = COMBINED_TABLES[fname]
    header_t, _ = read_tsv(path_t, [])
    header_r, _ = read_tsv(path_r, [])
    samples_t = combined_samples(header_t)
    samples_r = combined_samples(header_r)
    value_cols_t = [c for cols in samples_t.values() for c in cols.values()]
    value_cols_r = [c for cols in samples_r.values() for c in cols.values()]
    a = load_keyed(path_t, key_cols, value_cols_t, args.mod_decimals, mod_cols)
    b = load_keyed(path_r, key_cols, value_cols_r, args.mod_decimals, mod_cols)
    result = {
        "count": (len(a["keys"]), len(b["keys"])),
        "overlap": set_overlap(set(a["keys"]), set(b["keys"])),
        "id_kind": kind,
        "quant": {},
        "absent_quant_types": sorted({s.strip() for s in COMBINED_QUANT_SUFFIXES} - set(samples_t) - set(samples_r)),
    }
    one_sided = sorted(set(samples_t) ^ set(samples_r))
    if one_sided:
        raise InputError(f"quant type(s) {one_sided} present in only one of {path_t} / {path_r}")
    shared = a["keys"].keys() & b["keys"].keys()
    for qtype in samples_t:
        pooled: list[tuple[float | None, float | None]] = []
        per_exp = {}
        for exp_t, col_t in samples_t[qtype].items():
            exp_r = exp_map.get(exp_t, exp_t)
            if exp_r not in samples_r[qtype]:
                raise InputError(
                    f"{path_t}: sample {exp_t!r} ({qtype}) has no counterpart {exp_r!r} in {path_r}; "
                    f"reference samples: {sorted(samples_r[qtype])} (use --experiment-map)"
                )
            col_r = samples_r[qtype][exp_r]
            pairs = [(a["keys"][k][col_t], b["keys"][k][col_r]) for k in shared]
            pooled.extend(pairs)
            per_exp[exp_t] = quant_agreement(pairs, args.min_quant_n, log=True)
        result["quant"][qtype] = {"pooled": quant_agreement(pooled, args.min_quant_n, log=True), "per_experiment": per_exp}
    return result


def compare_tmt(path_t: Path, path_r: Path, args) -> dict:
    header_t, _ = read_tsv(path_t, [TMT_LAST_META_COLUMN])
    header_r, _ = read_tsv(path_r, [TMT_LAST_META_COLUMN])
    key_t, key_r = header_t[0], header_r[0]
    if key_t != key_r:
        raise InputError(f"first (key) column differs: {path_t} has {key_t!r}, {path_r} has {key_r!r}")
    ch_t = header_t[header_t.index(TMT_LAST_META_COLUMN) + 1 :]
    ch_r = header_r[header_r.index(TMT_LAST_META_COLUMN) + 1 :]
    if ch_t != ch_r:
        raise InputError(f"TMT channel columns differ: {path_t} {ch_t} vs {path_r} {ch_r}")
    a = load_keyed(path_t, [key_t], ch_t, args.mod_decimals, None)
    b = load_keyed(path_r, [key_r], ch_r, args.mod_decimals, None)
    # TMT-Integrator abundance/ratio values are already log2.
    return {
        "count": (len(a["keys"]), len(b["keys"])),
        "overlap": set_overlap(set(a["keys"]), set(b["keys"])),
        "channels": compare_values(a, b, ch_t, args.min_quant_n, log=False),
    }


def parse_params(path: Path) -> dict[str, str]:
    out = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if path.name == "fragger.params":
                line = line.split("#", 1)[0]
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sep = "=" if "=" in line else None
            if sep is None:
                continue
            key, value = line.split(sep, 1)
            value = value.strip().replace("\\:", ":").replace("\\\\", "\\")
            out[key.strip()] = value
    return out


def diff_params(path_t: Path, path_r: Path) -> dict:
    a, b = parse_params(path_t), parse_params(path_r)
    differing = {k: {"test": a.get(k), "reference": b.get(k)} for k in sorted(a.keys() | b.keys()) if a.get(k) != b.get(k)}
    return {
        "test_file": str(path_t),
        "reference_file": str(path_r),
        "differing": {k: v for k, v in differing.items() if not ENV_PARAM_PATTERN.search(k)},
        "differing_environment": {k: v for k, v in differing.items() if ENV_PARAM_PATTERN.search(k)},
    }


# --------------------------------------------------------------------------- #
# Gating
# --------------------------------------------------------------------------- #


def deviation(test: int, ref: int) -> float:
    return abs(test - ref) / max(test, ref, 1)


def evaluate(report: dict, args) -> list[dict]:
    checks: list[dict] = []

    def add(name: str, value: float | None, op: str, limit: float) -> None:
        if value is None:
            checks.append({"check": name, "value": None, "limit": limit, "status": "SKIP"})
            return
        ok = value <= limit if op == "<=" else value >= limit
        checks.append({"check": name, "value": value, "limit": limit, "op": op, "status": "PASS" if ok else "FAIL"})

    def add_spearman(name: str, metrics: dict) -> None:
        if metrics["spearman"] is None and (metrics["only_test"] or metrics["only_reference"]) and metrics["n"] < args.min_quant_n:
            # One side has values but too few overlap: a real disagreement.
            checks.append({"check": name, "value": None, "limit": args.min_spearman, "status": "FAIL", "note": metrics["note"]})
            return
        add(name, metrics["spearman"], ">=", args.min_spearman)

    for exp, res in report["experiments"].items():
        for what, (t, r) in res["counts"].items():
            add(f"{exp}: {what} count deviation", deviation(t, r), "<=", args.max_count_deviation)
        for what, ov in res["overlap"].items():
            add(f"{exp}: {what} Jaccard", ov["jaccard"], ">=", args.min_jaccard)
        add(f"{exp}: psm:spectra Jaccard", res["psm"]["spectra"]["jaccard"], ">=", args.min_jaccard)
        add(f"{exp}: shared-PSM same peptide (fraction)",
            None if res["psm"]["pct_same_peptide"] is None else res["psm"]["pct_same_peptide"] / 100, ">=", args.min_jaccard)
        for what, m in res["quant"].items():
            add_spearman(f"{exp}: {what} Spearman", m)
    for fname, res in report["combined"].items():
        add(f"{fname}: count deviation", deviation(*res["count"]), "<=", args.max_count_deviation)
        add(f"{fname}: {res['id_kind']} Jaccard", res["overlap"]["jaccard"], ">=", args.min_jaccard)
        for qtype, q in res["quant"].items():
            add_spearman(f"{fname}: {qtype} pooled Spearman", q["pooled"])
            for exp, m in q["per_experiment"].items():
                add_spearman(f"{fname}: {qtype} {exp} Spearman", m)
    for fname, res in report["tmt"].items():
        add(f"{fname}: row count deviation", deviation(*res["count"]), "<=", args.max_count_deviation)
        add(f"{fname}: row Jaccard", res["overlap"]["jaccard"], ">=", args.min_jaccard)
        for ch, m in res["channels"].items():
            add_spearman(f"{fname}: channel {ch} Spearman", m)
    return checks


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #


def fmt(value) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def quant_row(name: str, m: dict) -> str:
    return (f"| {name} | {m['n']} | {m['only_test']} | {m['only_reference']} | {fmt(m['pearson'])} | "
            f"{fmt(m['spearman'])} | {fmt(m['median_log2_ratio'])} | {fmt(m['pct_within_2fold'])} |")


QUANT_HEADER = ["| Metric | n shared | only test | only ref | Pearson | Spearman | median log2(test/ref) | % within 2-fold |",
                "| --- | --- | --- | --- | --- | --- | --- | --- |"]
OVERLAP_HEADER = ["| Set | test | reference | shared | only test | only ref | Jaccard |", "| --- | --- | --- | --- | --- | --- | --- |"]


def overlap_row(name: str, ov: dict) -> str:
    return f"| {name} | {ov['test']} | {ov['reference']} | {ov['shared']} | {ov['only_test']} | {ov['only_reference']} | {fmt(ov['jaccard'])} |"


def render_markdown(report: dict) -> str:
    lines = [
        "# FragPipe result comparison",
        "",
        f"**Verdict: {report['verdict']}**",
        "",
        f"- Test: `{report['test_dir']}`",
        f"- Reference: `{report['reference_dir']}`",
        f"- Thresholds: {', '.join(f'{k}={v}' for k, v in report['thresholds'].items())}",
        f"- Runtime: {report['runtime_seconds']:.1f} s",
        "",
    ]
    if report["unmatched_experiments"]["test"] or report["unmatched_experiments"]["reference"]:
        lines += ["## Unmatched experiments", "",
                  f"- Only in test: {report['unmatched_experiments']['test']}",
                  f"- Only in reference: {report['unmatched_experiments']['reference']}", ""]
    failed = [c for c in report["checks"] if c["status"] == "FAIL"]
    lines += ["## Failed checks", ""]
    if failed:
        lines += ["| Check | Value | Limit |", "| --- | --- | --- |"]
        lines += [f"| {c['check']} | {fmt(c['value'])} | {c.get('op', '')} {c['limit']} |" for c in failed]
    else:
        lines.append("None.")
    skipped = [c for c in report["checks"] if c["status"] == "SKIP"]
    lines += ["", f"Checks: {len(report['checks'])} total, {len(failed)} failed, {len(skipped)} skipped (metric undefined).", ""]

    if report["combined"]:
        lines += ["## Combined tables", ""] + OVERLAP_HEADER
        lines += [overlap_row(f"{f} ({r['id_kind']})", r["overlap"]) for f, r in report["combined"].items()]
        lines += [""] + QUANT_HEADER
        for f, r in report["combined"].items():
            for qtype, q in r["quant"].items():
                lines.append(quant_row(f"{f} {qtype} (pooled)", q["pooled"]))
        for f, r in report["combined"].items():
            if r["absent_quant_types"]:
                lines.append(f"\n{f}: quant type(s) absent on both sides: {r['absent_quant_types']}")
        lines.append("")
    for exp, r in report["experiments"].items():
        lines += [f"## Experiment {exp}", "", "| Count | test | reference | deviation |", "| --- | --- | --- | --- |"]
        lines += [f"| {k} | {t} | {ref} | {fmt(deviation(t, ref))} |" for k, (t, ref) in r["counts"].items()]
        lines += [""] + OVERLAP_HEADER + [overlap_row(k, ov) for k, ov in r["overlap"].items()]
        lines.append(overlap_row("psm: spectra", r["psm"]["spectra"]))
        lines += ["", f"Shared PSMs with same peptide: {fmt(r['psm']['pct_same_peptide'])} %; "
                  f"same modified peptide: {fmt(r['psm']['pct_same_modified_peptide'])} %.", ""]
        lines += QUANT_HEADER + [quant_row(f"score {k} (same-peptide PSMs, linear)", m) for k, m in r["psm"]["scores"].items()]
        lines += [quant_row(k, m) for k, m in r["quant"].items()]
        lines.append("")
    if report["tmt"]:
        lines += ["## TMT-Integrator tables", ""] + OVERLAP_HEADER
        lines += [overlap_row(f, r["overlap"]) for f, r in report["tmt"].items()]
        lines += ["", "Channel values are TMT-Integrator log2 values, compared without further transformation.", ""] + QUANT_HEADER
        for f, r in report["tmt"].items():
            lines += [quant_row(f"{f} {ch}", m) for ch, m in r["channels"].items()]
        lines.append("")
    if report["params"]:
        lines += ["## Parameter differences (informational, not gated)", ""]
        for fname, d in report["params"].items():
            lines += [f"### {fname}", ""]
            if d["differing"]:
                lines += ["| Key | test | reference |", "| --- | --- | --- |"]
                lines += [f"| {k} | `{v['test']}` | `{v['reference']}` |" for k, v in d["differing"].items()]
            else:
                lines.append("No non-environment differences.")
            lines += ["", f"Environment-specific keys differing (paths/threads): {sorted(d['differing_environment'])}", ""]
    if report["missing_params"]:
        lines += [f"Parameter files present on only one side (not diffed): {report['missing_params']}", ""]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


def parse_map(pairs: list[str]) -> dict[str, str]:
    out = {}
    for item in pairs:
        if "=" not in item:
            raise InputError(f"--experiment-map entry {item!r} is not TEST=REFERENCE")
        k, v = item.split("=", 1)
        out[k] = v
    return out


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compare msproteomics (test) FragPipe-style outputs against a FragPipe (reference) result directory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--test", required=True, type=Path, help="msproteomics result directory (searched recursively)")
    p.add_argument("--reference", required=True, type=Path, help="FragPipe result directory (searched recursively)")
    p.add_argument("--outdir", required=True, type=Path, help="directory for comparison.md and comparison.json")
    p.add_argument("--experiment-map", nargs="*", default=[], metavar="TEST=REF",
                   help="map test experiment names to reference names (default: identical names)")
    p.add_argument("--exclude", default=None, help="regex; files whose path matches are ignored during discovery")
    p.add_argument("--max-count-deviation", type=float, default=0.01, help="max |test-ref|/max(test,ref) for counts")
    p.add_argument("--min-jaccard", type=float, default=0.98,
                   help="min Jaccard for ID sets (also min fraction of shared PSMs with the same peptide)")
    p.add_argument("--min-spearman", type=float, default=0.99, help="min Spearman for intensities / TMT channels")
    p.add_argument("--min-quant-n", type=int, default=10, help="min shared quantified entries to compute correlations")
    p.add_argument("--mod-decimals", type=int, default=2, help="decimals used when matching modification masses")
    return p


def run(args) -> int:
    start = time.perf_counter()
    exclude = re.compile(args.exclude) if args.exclude else None
    found_t = discover(args.test, exclude)
    found_r = discover(args.reference, exclude)
    exp_map = parse_map(args.experiment_map)

    mapped = {t: exp_map.get(t, t) for t in found_t["experiments"]}
    matched = {t: r for t, r in mapped.items() if r in found_r["experiments"]}
    if (found_t["experiments"] or found_r["experiments"]) and not matched:
        raise InputError(
            f"no experiment directories match by name: test {sorted(found_t['experiments'])} vs "
            f"reference {sorted(found_r['experiments'])} (use --experiment-map)"
        )
    one_sided = sorted(set(found_t["combined"]) ^ set(found_r["combined"])) + sorted(set(found_t["tmt"]) ^ set(found_r["tmt"]))
    if one_sided:
        raise InputError(
            f"table(s) {one_sided} found on only one side "
            f"(test: {sorted(found_t['combined']) + sorted(found_t['tmt'])}, "
            f"reference: {sorted(found_r['combined']) + sorted(found_r['tmt'])})"
        )

    report = {
        "test_dir": str(args.test),
        "reference_dir": str(args.reference),
        "thresholds": {"max_count_deviation": args.max_count_deviation, "min_jaccard": args.min_jaccard,
                       "min_spearman": args.min_spearman, "min_quant_n": args.min_quant_n},
        "unmatched_experiments": {"test": sorted(set(found_t["experiments"]) - set(matched)),
                                  "reference": sorted(set(found_r["experiments"]) - set(matched.values()))},
        "experiments": {t: compare_experiment(t, r, found_t["experiments"][t], found_r["experiments"][r], args)
                        for t, r in sorted(matched.items())},
        "combined": {f: compare_combined(f, found_t["combined"][f], found_r["combined"][f], exp_map, args)
                     for f in COMBINED_FILES if f in found_t["combined"]},
        "tmt": {f: compare_tmt(found_t["tmt"][f], found_r["tmt"][f], args) for f in sorted(found_t["tmt"])},
        "params": {f: diff_params(found_t["params"][f], found_r["params"][f])
                   for f in PARAM_FILES if f in found_t["params"] and f in found_r["params"]},
        "missing_params": sorted(set(found_t["params"]) ^ set(found_r["params"])),
    }
    report["checks"] = evaluate(report, args)
    failed = any(c["status"] == "FAIL" for c in report["checks"]) or bool(
        report["unmatched_experiments"]["test"] or report["unmatched_experiments"]["reference"])
    report["verdict"] = "FAIL" if failed else "PASS"
    report["runtime_seconds"] = time.perf_counter() - start

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "comparison.json").write_text(json.dumps(report, indent=2, default=list) + "\n")
    (args.outdir / "comparison.md").write_text(render_markdown(report))
    n_fail = sum(c["status"] == "FAIL" for c in report["checks"])
    print(f"{report['verdict']}: {len(report['checks'])} checks, {n_fail} failed, "
          f"{report['runtime_seconds']:.1f} s -> {args.outdir / 'comparison.md'}")
    return EXIT_FAIL if failed else EXIT_PASS


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except InputError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":
    sys.exit(main())
