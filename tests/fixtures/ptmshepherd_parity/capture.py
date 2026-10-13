#!/usr/bin/env python3
"""PTM-Shepherd / PairScans parity fixtures from FragPipe 24.0 headless dry-run plans.

cases.json names each case: a workflow of tests/fixtures/fragpipe_workflows ("base"), optionally with
keys removed ("drop", full-key regexes) or set ("set"). FragPipe's `--dry-run --export-plan` writes every
tool's argv and config bytes without running anything (FragPipe branch export-plan-upstream-min).

    capture.py workflows <out dir>                 write the case workflows
    capture.sh <out dir> <plans dir>               (licensed image) export one plan per workflow
    capture.py expected <plans dir> <tools dir>    write expected/<case>.json from the plans

expected/<case>.json holds FragPipe's shepherd.config property block (the lines after the database /
dataset header) with the run's thread count and tools dir as the parser's @THREADS@ / @TOOLS_DIR@
tokens, the PairScans arguments after `<LC-MS file> <threads>` (CmdPairScans.java:121-126), O-Pair's
arguments after CMD.dll (CmdOPair.java:177-209) and MBG's after com.mbg.MBG (CmdMBGMatch.java:150-190)
with their run paths as the parser's tokens, and, when FragPipe refuses the workflow, its error message instead.
"""

import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKFLOWS = HERE.parent / "fragpipe_workflows"
THREADS = "4"  # capture.sh runs FragPipe with --threads 4


def cases() -> dict:
    return json.loads((HERE / "cases.json").read_text())


def case_workflow(case: dict) -> str:
    """The case's workflow text: the base workflow with `drop` keys removed and `set` keys replaced or appended."""
    values = dict(case.get("set", {}))
    lines = []
    for line in (WORKFLOWS / case["base"]).read_text().splitlines():
        key = line.split("=", 1)[0]
        if any(re.fullmatch(d, key) for d in case.get("drop", [])):
            continue
        if key in values:
            line = f"{key}={values.pop(key)}"
        lines.append(line)
    lines += [f"{k}={v}" for k, v in values.items()]
    return "\n".join(lines) + "\n"


def same(previous, tail, plan, tool):
    """One argv for every group: the group paths are tokens."""
    if previous not in (None, tail):
        raise ValueError(f"{plan}: {tool} arguments differ between groups")
    return tail


def expected(plan: Path, log: Path, tools_dir: str) -> dict:
    wd = re.escape(f"/tmp/wd/{plan.name.removesuffix('.plan.xml')}")  # capture.sh --workdir
    if not plan.exists():
        errors = [line.split(" ERROR - ", 1)[1] for line in log.read_text().splitlines() if " ERROR - " in line]
        return {"refusal": errors[0]}
    out = {"shepherd_config": None, "pairscans": None, "opair": None, "mbg": None}
    for step in ET.parse(plan).getroot().findall("step"):
        argv = [a.text or "" for a in step.find("argv")]
        if step.get("class") == "CmdPtmshepherd":
            config = next(c.text for c in step.findall("config") if c.get("path").endswith("/shepherd.config"))
            block = config.split("\n\n", 1)[1].replace(tools_dir, "@TOOLS_DIR@")
            block = re.sub(rf"^threads = {THREADS}$", "threads = @THREADS@", block, flags=re.M)
            out["shepherd_config"] = block.splitlines()
        elif step.get("class") == "CmdPairScans":
            tail = argv[argv.index("org.nesvilab.fragpipe.util.PairScans") + 3:]
            if out["pairscans"] not in (None, tail):
                raise ValueError(f"{plan}: PairScans arguments differ between LC-MS files")
            out["pairscans"] = tail
        elif step.get("class") == "CmdOPair":
            tail = [a.replace(tools_dir, "@TOOLS_DIR@") for a in argv[argv.index(next(a for a in argv if a.endswith("/CMD.dll"))) + 1:]]
            subs = [(rf"^-t {THREADS}$", "-t @THREADS@"), (rf"^-r {wd}/[^/]+/filelist_opair\.txt$", "-r @FILELIST@"),
                    (rf"^-r {wd}/filelist_opair\.txt$", "-r @FILELIST@"), (rf"^-s {wd}(/[^/]+)?/psm\.tsv$", "-s @PSM@"),
                    (rf"^-o {wd}(/[^/]+)?$", "-o @OUTDIR@")]
            for pattern, token in subs:
                tail = [re.sub(pattern, token, a) for a in tail]
            out["opair"] = same(out["opair"], tail, plan, "O-Pair")
        elif step.get("class") == "CmdMBGMatch":
            tail = [a.replace(tools_dir, "@TOOLS_DIR@") for a in argv[argv.index("com.mbg.MBG") + 1:]]
            tail = [re.sub(rf"^{wd}/[^/]+/psm\.tsv$", "@PSM@", a) for a in tail]
            tail = [re.sub(rf"^{wd}/fragpipe-files\.fp-manifest$", "@MANIFEST@", a) for a in tail]
            tail[tail.index("--numthreads") + 1] = "@THREADS@"
            tail[tail.index("--nopasef") + 1] = "@NOPASEF@"  # from the LC-MS files, filled by the MBG module
            out["mbg"] = same(out["mbg"], tail, plan, "MBG")
    text = json.dumps(out)
    if "/fragpipe_bin" in text or "/tmp/" in text:
        raise ValueError(f"{plan}: a run path is left in the fixture")
    return out


def main() -> None:
    command = sys.argv[1]
    if command == "workflows":
        out = Path(sys.argv[2])
        out.mkdir(parents=True, exist_ok=True)
        for name, case in cases().items():
            (out / f"{name}.workflow").write_text(case_workflow(case))
    elif command == "expected":
        plans, tools_dir = Path(sys.argv[2]), sys.argv[3]
        (HERE / "expected").mkdir(exist_ok=True)
        for name in cases():
            result = expected(plans / f"{name}.plan.xml", plans / f"{name}.log", tools_dir)
            (HERE / "expected" / f"{name}.json").write_text(json.dumps(result, indent=4) + "\n")
    else:
        raise SystemExit(f"unknown command {command}")


if __name__ == "__main__":
    main()
