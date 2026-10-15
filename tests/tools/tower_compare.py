#!/usr/bin/env python3
"""Pair a FragPipe results dir with a downloaded msproteomics (Tower/Fusion) outdir and run golden_compare.py pairs.

usage: tower_compare.py <golden files/work> <actual outdir> <json out> [--profile dia] [--only REGEX]
       tower_compare.py --ref-headless <ref outdir> <actual outdir> <json out> [--profile dia] [--only REGEX]

--ref-headless compares a FragPipe headless Tower run (<ref outdir>/fragpipe/results) without a golden manifest.json:
the roots are read from its fragpipe.job (the /tmp/nxf.<id> task dir, results/ = workdir, raw_files/ = inputs).

<golden files/work> sits two levels below a manifest.json whose "roots" name the recorded run's directories
({WORKDIR}, {LCMS_DIR}, {DB_DIR}, {TMP}); they, and every Nextflow task dir the actual files name, map to <RUN>.

Pairing: same basename. When the basename occurs more than once on the actual side, the file whose parent dir name
equals the golden file's parent dir name wins; then the stage that rewrites the final tables (STAGE_ORDER) wins.
Golden files with no actual counterpart are listed as MISSING (not compared).
"""
import argparse
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# A Nextflow task dir as FragPipe-in-Nextflow files record it: <fusion root>/<any work dir>/<xx>/<30 hex>, or the same
# dir as an S3 URL (a path Groovy writes, e.g. the run-level fragger*.params database_name)
TASK_DIR = rb"(?:/fusion/s3/|s3://)[a-z0-9][a-z0-9.-]*/(?:[^\s\t\"'<>/]+/)*?[0-9a-f]{2}/[0-9a-f]{30}"
# METAPROTEOMICS runs FP-Meta in <task dir>/project (project_resolved before the copy was dropped), the counterpart
# of FragPipe's work dir
PROJECT_RESOLVED = TASK_DIR + rb"/project(?:_resolved)?(?![A-Za-z0-9_.-])"
# scratch = true task: the task runs in /tmp/nxf.<id>; FragPipe headless: results/ = workdir, raw_files/ = inputs
NXF_TMP = rb"/tmp/nxf\.[A-Za-z0-9]+(?:/results|/raw_files)?"
# Input dirs named as S3 URLs (the "/" after "s3:/") or Fusion paths: <bucket>/.../inputs[/raw|/d]/, or in-place spectra
_BUCKET = rb"(?:/fusion/s3|(?<=s3:/))/[a-z0-9][a-z0-9.-]*"
INPUT_DIRS = (
    _BUCKET + rb"/[^\s\t\"'<>]*?/inputs(?:/raw|/d)?(?=/)"
    rb"|" + _BUCKET + rb"/[^\s\t\"'<>]*/(?=[^/\s]+\.(?:raw|d)\b)"
)
PATTERN = re.compile(PROJECT_RESOLVED + b"|" + TASK_DIR + b"|" + NXF_TMP + b"|" + INPUT_DIRS)
TASK_DIR_RE = re.compile(TASK_DIR)

# Stages that rewrite psm/ion/peptide/protein.tsv after Philosopher, earliest preference first
STAGE_ORDER = ("ionquant_isobaric", "ionquant_ms1", "ionquant", "mbg", "opair", "ptmshepherd")


def is_task_dir(prefix: str) -> bool:
    """True for a task dir prefix (Philosopher's --temp lives below it in phi_temp/)."""
    return TASK_DIR_RE.fullmatch(prefix.encode()) is not None


def choose_candidate(cands: list[Path], golden: Path) -> list[Path]:
    """Narrow the same-basename files of the actual side to the one that pairs with `golden` (one item if resolved)."""
    if len(cands) > 1:
        cands = [c for c in cands if c.parent.name == golden.parent.name] or cands
    if len(cands) > 1 and re.fullmatch(r"fragger.*\.params", golden.name):
        # FragPipe's run-level params (CmdMsfragger.java:444-447) pair with msfragger/<name>, not the per-file copies
        cands = [c for c in cands if c.parent.name == "msfragger"] or cands
    for stage in STAGE_ORDER:
        if len(cands) > 1:
            cands = [c for c in cands if stage in c.parts] or cands
    return cands


def headless_roots(results: Path) -> dict:
    """manifest roots of a FragPipe headless task, from the /tmp/nxf.<id> dir its fragpipe.job names."""
    m = re.search(rb"/tmp/nxf\.[A-Za-z0-9]+", (results / "fragpipe.job").read_bytes())
    if m is None:
        raise ValueError(f"{results / 'fragpipe.job'}: no /tmp/nxf.<id> task dir named")
    nxf = m.group().decode()
    return {"{WORKDIR}": nxf + "/results", "{LCMS_DIR}": nxf + "/raw_files", "{DB_DIR}": nxf,
            "{FRAGPIPE}": "/fragpipe_bin/fragpipe-24.0/fragpipe-24.0", "{TMP}": "/tmp"}


def pair_files(g_root: Path, a_root: Path, only: str = "."):
    by_name = defaultdict(list)
    for p in a_root.rglob("*"):
        if p.is_file() and "pipeline_info" not in p.parts and "multiqc" not in p.parts:
            by_name[p.name].append(p)
    pairs, missing = [], []
    for g in sorted(p for p in g_root.rglob("*") if p.is_file()):
        rel = g.relative_to(g_root).as_posix()
        if not re.search(only, rel):
            continue
        cands = choose_candidate(by_name.get(g.name, []), g)
        if len(cands) == 1:
            pairs.append((rel, g, cands[0]))
        else:
            missing.append((rel, len(cands)))
    return pairs, missing


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-headless", metavar="REF_OUTDIR")
    ap.add_argument("paths", nargs="+", help="<golden files/work> <actual outdir> <json out>, or with --ref-headless <actual outdir> <json out>")
    ap.add_argument("--profile", default="dda")
    ap.add_argument("--only", default=r".")
    a = ap.parse_args(argv)
    if len(a.paths) != (2 if a.ref_headless else 3):
        ap.error("expected <golden files/work> <actual outdir> <json out>, or with --ref-headless <actual outdir> <json out>")
    if a.ref_headless:
        g_root = Path(a.ref_headless) / "fragpipe/results"
        a_root, a.json_out = Path(a.paths[0]), a.paths[1]
    else:
        g_root, a_root, a.json_out = Path(a.paths[0]), Path(a.paths[1]), a.paths[2]
    pairs, missing = pair_files(g_root, a_root, a.only)
    manifest = {"roots": headless_roots(g_root)} if a.ref_headless else json.loads((g_root.parent.parent / "manifest.json").read_text())
    exp_prefixes = sorted({v for k, v in manifest["roots"].items() if k != "{FRAGPIPE}"}, key=len, reverse=True)
    act_prefixes = set()
    for _, _, c in pairs:
        with open(c, "rb") as fh:
            act_prefixes.update(m.decode() for m in PATTERN.findall(fh.read()))
    gc = Path(__file__).resolve().parent / "golden_compare.py"
    args = [sys.executable, str(gc), "pairs", "--profile", a.profile, "--json", a.json_out]
    for x in exp_prefixes:
        args += ["--expected-prefix", x]
    for x in sorted(act_prefixes, key=len, reverse=True):
        args += ["--actual-prefix", x.rstrip("/")]
    for x in sorted(act_prefixes):
        if is_task_dir(x):  # philosopher --temp of the task (phi_temp/<uuid>/)
            args += ["--actual-tmp", x + "/phi_temp"]
    if any(x.startswith("/tmp/nxf.") for x in act_prefixes):  # scratch = true task: Philosopher --temp is /tmp/<uuid>
        args += ["--actual-tmp", "/tmp"]
    if "{TMP}" in manifest["roots"]:
        args += ["--expected-tmp", manifest["roots"]["{TMP}"]]
    for rel, g, c in pairs:
        args += ["--pair", rel, str(g), str(c)]
    r = subprocess.run(args, capture_output=True, text=True)
    print(r.stdout[-200:] if r.returncode not in (0, 1) else "", r.stderr[-3000:], file=sys.stderr)
    res = json.loads(Path(a.json_out).read_text())
    res["missing_on_actual_side"] = [{"path": m, "candidates": n} for m, n in missing]
    Path(a.json_out).write_text(json.dumps(res, indent=1) + "\n")
    for x in res["results"]:
        print(f"{x['status']:6s} {x['path']}  {x.get('rule', '')}  {x.get('detail', '')[:160]}")
    print("counts:", res["counts"], "missing:", len(missing))
    for m, n in missing:
        print(f"MISSING {m} ({n} candidates)")


if __name__ == "__main__":
    main()
