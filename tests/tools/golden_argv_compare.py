#!/usr/bin/env python3
"""Compare one logged module command with a recorded FragPipe command (no manifest step).

For commands FragPipe runs inside a script (e.g. ``msfragger_pep_split.py``), the manifest holds only the
script's argv. The inner commands are recorded once from the FragPipe script itself
(``tests/fixtures/golden/*.json``: ``{"stem": ..., "commands": {key: [argv...]}}``). Both sides go through
the accepted normalisations of ``golden_step_check`` (``-Duser.home``, ``-Xmx<N>``, licence path) plus the
location tokens that differ between FragPipe's run and a Nextflow task directory:

- paths under the FragPipe install become ``{FRAGPIPE}``
- the run stem becomes ``{STEM}`` (the recorded run and the compared run name their files differently)
- a path ending in ``fragger.params`` (or FragPipe's per data type ``fragger_{dda,dia,gpfdia,dda_plus}.params``,
  CmdMsfragger.java:444-447) becomes ``{PARAMS}`` (FragPipe's tempdir copy, the module's staged copy)
- a path ending in ``.mzML`` or ``.mzBIN_calibrated`` becomes ``{SPECTRA:<file name>}``

The last two drop the directory, so each side names the directories its files may be in: the recorded
fixture's ``dirs`` and ``--actual-dir`` (the module task's own directories). A file in any other directory
keeps its full path and so is a DIFF; a bare file name is in the command's working directory and allowed.

Exit code 0 only if the normalised argv are equal.
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
import golden_step_check as gsc  # noqa: E402


def normalise(argv: list[str], fragpipe: str, stem: str, dirs: list[str]) -> list[str]:
    """Accepted normalisations, then the location tokens, on one argv; ``dirs`` are the allowed directories."""
    allowed = {d.rstrip("/").replace(stem, "{STEM}") for d in dirs} | {""}
    out = []
    for tok in gsc.accepted_normalisations(argv):
        tok = tok.replace(stem, "{STEM}")
        directory, name = os.path.split(tok)
        if tok.startswith(fragpipe.rstrip("/") + "/"):
            tok = "{FRAGPIPE}/" + tok[len(fragpipe.rstrip("/")) + 1 :]
        elif directory.rstrip("/") not in allowed:
            pass  # a file outside the side's own directories: compared with its full path
        elif re.fullmatch(r"(.*_)?fragger(_dda|_dia|_gpfdia|_dda_plus)?\.params", name):
            tok = "{PARAMS}"
        elif name.endswith((".mzML", ".mzBIN_calibrated")):
            tok = "{SPECTRA:" + name + "}"
        out.append(tok)
    return out


def compare(expected: list[str], actual: list[str], fragpipe: str, expected_stem: str, actual_stem: str, expected_dirs: list[str], actual_dirs: list[str]) -> dict:
    exp = normalise(expected, fragpipe, expected_stem, expected_dirs)
    act = normalise(actual, fragpipe, actual_stem, actual_dirs)
    if exp == act:
        return {"status": "MATCH", "expected": exp, "actual": act}
    return {"status": "DIFF", "expected": exp, "actual": act, "token_diff": [d for d in difflib.ndiff(exp, act) if d[:1] in "+-"]}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--expected", type=Path, required=True, help="recorded commands JSON")
    ap.add_argument("--key", required=True, help="command name in the recorded JSON")
    ap.add_argument("--actual-cmd", type=Path, required=True, help="file whose first line holds the module's command")
    ap.add_argument("--actual-stem", required=True, help="run stem of the compared module run")
    ap.add_argument("--actual-fragpipe", help="FragPipe install root the module used (default: taken from the command, golden_step_check.fragpipe_root)")
    ap.add_argument("--actual-dir", action="append", default=[], help="a directory of the module task the command may name files in (repeatable)")
    a = ap.parse_args(argv)
    fx = json.loads(a.expected.read_text(encoding="utf-8"))
    line = a.actual_cmd.read_text(encoding="utf-8").splitlines()[0]
    actual = shlex.split(line)
    res = compare(fx["commands"][a.key], actual, a.actual_fragpipe or gsc.fragpipe_root(line), fx["stem"], a.actual_stem, fx["dirs"], a.actual_dir)
    print(f"ARGV {res['status']} {a.key}")
    if res["status"] != "MATCH":
        print("  expected: " + " ".join(res["expected"]))
        print("  actual:   " + " ".join(res["actual"]))
        for d in res["token_diff"]:
            print(f"  {d}")
    return 0 if res["status"] == "MATCH" else 1


if __name__ == "__main__":
    sys.exit(main())
