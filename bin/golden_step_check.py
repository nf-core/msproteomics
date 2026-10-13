#!/usr/bin/env python3
"""Check one msproteomics module run against FragPipe golden steps.

Given a golden set (``bin/golden_capture.py``) and the manifest indices of the
FragPipe steps the module reimplements, this checks:

(a) argv: the module's command line for each step equals FragPipe's argv after
    path normalization. Both sides become step-relative: the step's work dir is
    ``{CWD}``, the FragPipe install is ``{FRAGPIPE}``, external inputs are
    ``{INPUT:<name>}``. Nextflow stages every input into the task dir, so a
    module token naming a staged input (``{CWD}/<name>`` or bare ``<name>``)
    takes the location token FragPipe used for that input. Nothing else is
    normalized: argument order, extra or missing flags and different
    executables are reported as DIFF.

(b) outputs: every output FragPipe left on disk for these steps is compared
    with the module's file of the same name using ``bin/golden_compare.py``
    rules. Outputs FragPipe deleted, moved or overwrote later are reported as
    UNVERIFIABLE (no golden copy exists); a step with no MATCHed output at
    all fails as NO_VERIFIED_OUTPUT. Module files in ``--actual-outdir``
    that match no step output are EXTRA, unless listed with ``--ignore-actual``.

Exit code 0 only if every argv and every verifiable file MATCHes.
"""

from __future__ import annotations

import argparse
import difflib
import json
import os
import re
import shlex
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import golden_compare as gc  # noqa: E402

CWD = "{CWD}"
INPUT_TOKEN_RE = re.compile(r"\{INPUT:([^}]+)\}")


def normpaths(tok: str) -> str:
    """Resolve ``.``/``..`` in each absolute path of a token (``:``-separated lists too)."""
    return ":".join(os.path.normpath(p) if p.startswith("/") else p for p in tok.split(":"))


def load_golden(golden: Path) -> dict:
    m = json.loads((golden / "manifest.json").read_text(encoding="utf-8"))
    if "golden" not in m:
        raise ValueError(f"{golden}/manifest.json has no 'golden' block; capture it with bin/golden_capture.py")
    return m


def golden_path(golden: Path, manifest: dict, placeholder_path: str) -> Path:
    for key, rel in manifest["golden"]["layout"].items():
        if placeholder_path == key or placeholder_path.startswith(key + "/"):
            return golden / rel / placeholder_path[len(key) + 1 :]
    raise ValueError(f"{placeholder_path} is outside the captured layout")


def _step_relative(tok: str, cwd: str | None) -> str:
    if cwd is not None and (tok == cwd or tok.startswith(cwd + "/")):
        return CWD + tok[len(cwd) :]
    return tok.replace(cwd + "/", CWD + "/") if cwd else tok


def expected_argv(step: dict) -> list[str]:
    """FragPipe argv in step-relative form; bare cwd file names become {CWD}/<name>."""
    cwd = step["cwd"]
    local = {
        os.path.basename(f["path"])
        for f in step["inputs"] + step["outputs"]
        if cwd is not None and os.path.dirname(f["path"]) == cwd
    }
    out = []
    for t in step["argv"]:
        t = _step_relative(t, cwd)
        out.append(f"{CWD}/{t}" if t in local else t)
    return out


FRAGPIPE_SUBDIR_RE = re.compile(r"/(?:tools|lib)/")


def fragpipe_root(cmdline: str) -> str:
    """FragPipe install root of a module command. The module gets its tool paths from
    fragpipeTools(task.ext.fragpipe_tools_dir): <root>/tools/... and FragPipe's own <root>/lib/*, where the
    root is a FragPipe install dir (named fragpipe*, e.g. /fragpipe_bin/fragpipe-24.0/fragpipe-24.0). Other
    tools/ or lib/ dirs (an absolute /usr/lib/jvm/.../java) are not taken. Every FragPipe path of the command
    must name the same root; none, or two different ones, is an error."""
    roots = []
    for t in shlex.split(cmdline, posix=True):
        for part in t.split(":"):
            for m in FRAGPIPE_SUBDIR_RE.finditer(part):
                root = part[:m.start()]
                if root.startswith("/") and os.path.basename(root).lower().startswith("fragpipe"):
                    roots.append(root)
                    break
    if not roots:
        raise ValueError(f"no <fragpipe*>/tools/ or <fragpipe*>/lib/ path in the module command: {cmdline[:200]}")
    if len(set(roots)) > 1:
        raise ValueError(f"module command names several FragPipe roots {sorted(set(roots))}: {cmdline[:200]}")
    return roots[0]


def actual_argv(cmdline: str, step: dict, actual_cwd: str, actual_fragpipe: str | None, actual_workdir: str | None = None, actual_lcms_dir: str | None = None, run_inputs: list[dict] | None = None, input_dir: str | None = None) -> list[str]:
    """Module argv in the same form. The logged command line is split with ``shlex``
    (POSIX shell rules), so quoted arguments and empty arguments (``""``) are kept,
    like the empty arguments ``split_command`` keeps on the FragPipe side."""
    by_name = {}
    # the run's own LC-MS files: a step whose inputs the manifest does not record (io_convention_not_implemented)
    # still names them as {INPUT:<name>}; a module token naming the staged file takes that location
    for f in run_inputs or []:
        m = INPUT_TOKEN_RE.fullmatch(f["path"])
        by_name[m.group(1) if m else os.path.basename(f["path"])] = f["path"]
    for f in step["inputs"] + step["outputs"]:
        loc = _step_relative(f["path"], step["cwd"])
        m = INPUT_TOKEN_RE.fullmatch(f["path"])
        by_name[m.group(1) if m else os.path.basename(f["path"])] = loc
    pairs = [(actual_cwd.rstrip("/"), CWD), ((actual_fragpipe or fragpipe_root(cmdline)).rstrip("/"), "{FRAGPIPE}")]
    if actual_workdir:  # module dir that stands for FragPipe's work dir
        pairs.append((actual_workdir.rstrip("/"), "{WORKDIR}"))
    if actual_lcms_dir:  # staged LC-MS directory that stands for FragPipe's input dir
        pairs.append((actual_lcms_dir.rstrip("/"), "{LCMS_DIR}"))
    if input_dir:  # module dir of the staged run inputs (LC-MS files, database); names in run_inputs take their {INPUT:..} token
        pairs.append((input_dir.rstrip("/"), "{INPUT_DIR}"))
    roots = sorted(pairs, key=lambda kv: len(kv[0]), reverse=True)
    out = []
    for t in shlex.split(cmdline, posix=True):
        t = normpaths(t)
        for root, token in roots:
            t = t.replace(root + "/", token + "/") if t != root else token
        name = None
        for root in (CWD, "{WORKDIR}", "{INPUT_DIR}"):
            if t.startswith(root + "/") and "/" not in t[len(root) + 1 :]:
                name = t[len(root) + 1 :]
        if name is None and "/" not in t:
            name = t
        if name is not None and name in by_name:
            t = by_name[name]
        out.append(t)
    return out


def rule_profile(manifest: dict) -> str:
    """Comparison rule profile from the run's data types: DIA rules only for DIA runs."""
    types = {f["data_type"] for f in manifest["lcms_files"]}
    if types & {"DIA", "DIA-Lib", "GPF-DIA", "DIA-Quant"}:
        return "dia"
    return "dda"


XMX_RE = re.compile(r"^-Xmx\d+[GgMm]$")
# FileMove's -cp is <fragpipe jar>:<the jar FragPipe's JVM loaded org.apache.commons.io.FileUtils from>
# (ToolingUtils.java:144-148, FragPipe 24.0); with lib/* on the classpath that is commons-io-2.18.0.jar or a
# batmass-io jar, depending on the install's lib/ contents and class-loading order, not on the step
FILEMOVE_CP_RE = re.compile(r"^(.+/lib/fragpipe-[^/:]+\.jar):(.+/lib/)(?:commons-io-2\.18\.0|batmass-io-[0-9.]+)\.jar$")
# Placeholder of FileMove's commons-io jar (a plain string: nf-core lint reads escaped f-string braces as Jinja)
COMMONS_IO_TOKEN = "{COMMONS_IO_JAR}"
FILEMOVE_CLASS = "org.nesvilab.utils.FileMove"
# MSFragger, IonQuant and diaTracer jars carry an edition word before the version in some FragPipe
# distributions; FragPipe finds either name (TabConfig.java:176-178, :188-209), so the word is not compared
EDITION_JAR_RE = re.compile(r"\b(MSFragger|IonQuant|diaTracer)-[A-Za-z][A-Za-z0-9]*-(\d+(?:\.\d+)+\.jar)")


def accepted_normalisations(argv: list[str]) -> list[str]:
    """Differences accepted between a module's argv and FragPipe's, applied to both sides.

    - ``-Duser.home=...`` and ``-Djava.io.tmpdir=...`` are dropped: the modules pass them per Java
      command (user decision; the container user has no passwd entry and the task dir is the tmp dir).
    - ``-Xmx<N>G`` becomes ``-Xmx<N>``: the heap comes from the task memory (the module subtracts
      headroom); evidence: the IonQuant and TMT-Integrator golden outputs match byte for byte with it.
    - The licence path after ``--license`` / in ``--license=`` becomes ``{LICENSE}``, the placeholder
      ``bin/golden_capture.py`` writes into the manifest (the real path is never stored).
    - In a FileMove command only, the second ``-cp`` entry (``lib/commons-io-2.18.0.jar`` or
      ``lib/batmass-io-<v>.jar``) becomes ``lib/{COMMONS_IO_JAR}`` (FILEMOVE_CP_RE).
    - ``MSFragger-<edition>-<v>.jar`` (also IonQuant, diaTracer) becomes ``MSFragger-<v>.jar`` (EDITION_JAR_RE).
    """
    out: list[str] = []
    argv = [EDITION_JAR_RE.sub(r"\1-\2", t) for t in argv]
    for i, t in enumerate(argv):
        if t.startswith(("-Duser.home=", "-Djava.io.tmpdir=")):
            continue
        m = FILEMOVE_CP_RE.match(t)
        if m and i > 0 and argv[i - 1] == "-cp" and i + 1 < len(argv) and argv[i + 1] == FILEMOVE_CLASS:
            t = f"{m.group(1)}:{m.group(2)}{COMMONS_IO_TOKEN}"
        elif XMX_RE.match(t):
            t = "-Xmx<N>"
        elif t.startswith("--license=") and t != "--license={LICENSE}":
            t = "--license={LICENSE}"
        elif i > 0 and argv[i - 1] == "--license" and t != "{LICENSE}":
            t = "{LICENSE}"
        out.append(t)
    return out


def check_argv(exp: list[str], act: list[str]) -> dict:
    exp, act = accepted_normalisations(exp), accepted_normalisations(act)
    if exp == act:
        return {"status": "MATCH", "expected": exp, "actual": act}
    diff = [d for d in difflib.ndiff(exp, act) if d[:1] in "+-"]
    return {"status": "DIFF", "expected": exp, "actual": act, "token_diff": diff}


def resolve_workspace_steps(steps: list[dict]) -> None:
    """A Philosopher step with no recorded output file (e.g. PhilosopherDbAnnotate) writes only the workspace
    (.meta) of its dir, which FragPipe's workspace clean deletes. The next checked step in the same dir reads
    that workspace, so a MATCHing next step verifies it (VIA_CONSUMER)."""
    for e in steps:
        if e["files"] or not e["cmd_name"].startswith("Philosopher"):
            continue
        later = [x for x in steps if x["step_index"] > e["step_index"] and x["cwd"] == e["cwd"]]
        nxt = min(later, key=lambda x: x["step_index"]) if later else None
        if nxt and nxt["argv"]["status"] == "MATCH" and any(x["status"] == "MATCH" for x in nxt["files"]) and all(x["status"] in ("MATCH", "UNVERIFIABLE", "VIA_CONSUMER") for x in nxt["files"]):
            e["files"].append({"path": f"{e['cwd']}/.meta", "status": "VIA_CONSUMER", "detail": f"workspace only; next step {nxt['step_index']} in the same dir MATCH"})


def actual_prefixes(a: argparse.Namespace, fp_dirs: list[str]) -> list[str]:
    """Run-specific path prefixes of the module side: its cwd, FragPipe's run dirs, the module work dir and the
    dir the module's database was read from (--actual-db-dir: FragPipe's {DB_DIR}, e.g. the golden db dir that a
    test points the input pepXML's search_database at)."""
    return [a.actual_cwd, *fp_dirs] + [d for d in (a.actual_workdir, getattr(a, "actual_db_dir", None)) if d]


def run(a: argparse.Namespace) -> dict:
    golden = a.golden
    m = load_golden(golden)
    steps = {s["step_index"]: s for s in m["steps"]}
    sel = [steps[i] for i in a.steps]
    roots = m["roots"]
    profile = a.profile or rule_profile(m)
    report_profile = profile
    exp_side = gc.Side(prefixes=[], tmp_roots=[roots["{TMP}"]] if "{TMP}" in roots else [])
    fp_dirs = [roots[k] for k in ("{WORKDIR}", "{LCMS_DIR}", "{DB_DIR}")]
    act_side = gc.Side(prefixes=actual_prefixes(a, fp_dirs), tmp_roots=a.actual_tmp)
    argv_cmds = dict(x.split("=", 1) for x in a.actual_argv)
    run_inputs = list(m["lcms_files"]) + ([{"path": m["database"]}] if m.get("database") else [])

    report: dict = {"golden": str(golden), "profile": report_profile, "steps": []}
    seen_actual: set[str] = set()
    for s in sel:
        cwd_abs = s["cwd"]
        for k, v in roots.items():
            if cwd_abs is not None:
                cwd_abs = cwd_abs.replace(k, v)
        exp_side.prefixes = [p for p in [cwd_abs, *fp_dirs] if p]
        entry: dict = {"step_index": s["step_index"], "name": s["name"], "cmd_name": s["cmd_name"], "cwd": s["cwd"]}
        key = str(s["step_index"])
        if key not in argv_cmds:
            raise ValueError(f"no --actual-argv given for step {key} ({s['name']})")
        cmdline = Path(argv_cmds[key]).read_text(encoding="utf-8").strip()
        entry["argv"] = check_argv(expected_argv(s), actual_argv(cmdline, s, a.actual_cwd, a.actual_fragpipe, a.actual_workdir, a.actual_lcms_dir, run_inputs, a.actual_input_dir))
        files = []
        for o in [] if a.argv_only else s["outputs"]:
            name = os.path.basename(o["path"])
            if o["state"] != "present":
                files.append({"path": o["path"], "status": "UNVERIFIABLE", "detail": f"FragPipe output {o['state']}; no golden copy"})
                seen_actual.add(name)
                continue
            pe = golden_path(golden, m, o["path"])
            pa = a.actual_outdir / name
            seen_actual.add(name)
            if not pa.is_file():
                files.append({"path": o["path"], "status": "MISSING", "detail": f"module did not write {name}"})
                continue
            rel = o["path"].split("/", 1)[1] if o["path"].startswith("{WORKDIR}/") else name
            rewritten = any(t["step_index"] > s["step_index"] and any(x["path"] == o["path"] for x in t["outputs"])
                            for t in m["steps"])
            r = gc.compare_files(pe, pa, rel, exp_side, act_side, profile, rewritten)
            files.append({"path": o["path"], "status": r.status, "rule": r.rule, "detail": r.detail})
        entry["files"] = files
        report["steps"].append(entry)
    # An output FragPipe overwrote later (e.g. PeptideProphet's interact.pep.xml, rewritten by RewritePepxml)
    # has no golden copy; it is verified through a checked step that reads it when that step's argv and
    # every file MATCH (VIA_CONSUMER). Without such a consumer it stays UNVERIFIABLE.
    by_index = {e["step_index"]: e for e in report["steps"]}
    for e in report["steps"]:
        for f in e["files"]:
            if f["status"] != "UNVERIFIABLE":
                continue
            consumers = [by_index[t["step_index"]] for t in sel if t["step_index"] in by_index and any(i["path"] == f["path"] for i in t["inputs"])]
            if consumers and all(c["argv"]["status"] == "MATCH" and any(x["status"] == "MATCH" for x in c["files"]) and all(x["status"] in ("MATCH", "UNVERIFIABLE", "VIA_CONSUMER") for x in c["files"]) for c in consumers):
                f["status"] = "VIA_CONSUMER"
                f["detail"] += "; consumer step(s) " + ", ".join(str(c["step_index"]) for c in consumers) + " MATCH"
    resolve_workspace_steps(report["steps"])
    for e in report["steps"]:
        n_ok = sum(1 for f in e["files"] if f["status"] in ("MATCH", "VIA_CONSUMER"))
        if n_ok == 0 and not a.argv_only:
            # A step whose outputs cannot be checked (none, or all UNVERIFIABLE) is not verified.
            e["files"].append({"path": "-", "status": "NO_VERIFIED_OUTPUT", "detail": "no output of this step could be compared with FragPipe"})
        e["verified_files"] = n_ok
    # Module files named like an input of a checked step are FragPipe intermediates
    # (e.g. Percolator *_psms.tsv read by PercolatorOutputToPepXML), not extras.
    step_inputs = {os.path.basename(i["path"]) for st in sel for i in st["inputs"]}
    extras = []
    for p in [] if a.argv_only else sorted(a.actual_outdir.iterdir()):
        if not p.is_file() or p.name in seen_actual:
            continue
        if p.name in a.ignore_actual:
            extras.append({"path": p.name, "status": "SKIP", "detail": "--ignore-actual"})
        elif p.name in step_inputs:
            extras.append({"path": p.name, "status": "SKIP", "detail": "input of a checked FragPipe step"})
        else:
            extras.append({"path": p.name, "status": "EXTRA"})
    report["actual_extra_files"] = extras
    bad = [e for e in report["steps"] if e["argv"]["status"] != "MATCH"]
    bad += [f for e in report["steps"] for f in e["files"] if f["status"] in ("DIFF", "MISSING", "NO_VERIFIED_OUTPUT")]
    bad += [x for x in extras if x["status"] == "EXTRA"]
    report["result"] = "FAIL" if bad else "PASS"
    return report


def print_report(rep: dict) -> None:
    for e in rep["steps"]:
        print(f"STEP {e['step_index']} {e['name']} [{e['cmd_name']}] cwd={e['cwd']}")
        print(f"  ARGV {e['argv']['status']}")
        if e["argv"]["status"] != "MATCH":
            print("    expected: " + " ".join(e["argv"]["expected"]))
            print("    actual:   " + " ".join(e["argv"]["actual"]))
            for d in e["argv"]["token_diff"]:
                print(f"    {d}")
        for f in e["files"]:
            print(f"  {f['status']}\t{f['path']}\t{f.get('rule', '')}\t{f.get('detail', '')}")
    for x in rep["actual_extra_files"]:
        print(f"  {x['status']}\t(actual) {x['path']}")
    print(f"PROFILE {rep['profile']}; verified files per step: " + ", ".join(f"{e['step_index']}={e['verified_files']}" for e in rep["steps"]))
    print(f"RESULT {rep['result']}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--golden", type=Path, required=True, help="golden/<dataset>/<workflow> dir")
    ap.add_argument("--steps", type=int, nargs="+", required=True, help="manifest step_index values covered by the module")
    ap.add_argument("--actual-cwd", required=True, help="dir the module ran its commands in (task work dir)")
    ap.add_argument("--actual-outdir", type=Path, required=True, help="dir holding the module's output files")
    ap.add_argument("--actual-fragpipe", help="FragPipe install root the module used (default: taken from each module command, fragpipe_root)")
    ap.add_argument("--actual-argv", action="append", default=[], metavar="STEP=FILE", help="file with the module's command line for a step")
    ap.add_argument("--actual-tmp", action="append", default=[], help="tmp root of the module run (for prot.xml)")
    ap.add_argument("--actual-workdir", help="module dir standing for FragPipe's work dir (optional)")
    ap.add_argument("--actual-lcms-dir", help="module dir of the staged LC-MS files, standing for FragPipe's {LCMS_DIR} (optional)")
    ap.add_argument("--actual-input-dir", help="module dir of the staged run inputs (LC-MS files, database), standing for FragPipe's {INPUT:<name>} locations (optional)")
    ap.add_argument("--actual-db-dir", help="dir the module's database was read from, standing for FragPipe's {DB_DIR} in outputs (optional)")
    ap.add_argument("--argv-only", action="store_true", help="check the argv only: no output files, no extra-file check (outputs are compared elsewhere, or the manifest records none: the effect lives in the Philosopher .meta workspace or stdout)")
    ap.add_argument("--ignore-actual", action="append", default=[], help="module-only diagnostic file name (e.g. percolator.log)")
    ap.add_argument("--profile", choices=gc.PROFILES, help="comparison rule profile (default: from the manifest's data types)")
    ap.add_argument("--json", type=Path, help="write the report as JSON")
    a = ap.parse_args(argv)
    rep = run(a)
    print_report(rep)
    if a.json:
        a.json.write_text(json.dumps(rep, indent=2) + "\n", encoding="utf-8")
    return 0 if rep["result"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
