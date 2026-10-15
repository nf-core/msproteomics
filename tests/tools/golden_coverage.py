#!/usr/bin/env python3
"""List every FragPipe step of a golden manifest and the golden nf-test covering it.

A golden test covers a step when its ``tests/golden.nf.test`` selects the step
from the manifest with ``it.cmd_name == '<cmd_name>'`` (optionally also
``it.cwd == '<cwd>'``). Prints a Markdown table to stdout.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

SELECT_RE = re.compile(r"it\.cmd_name == '([^']+)'(?: && it\.cwd == '([^']+)')?")


def find_tests(repo: Path) -> list[tuple[str, str, str | None]]:
    """(test path, cmd_name, cwd or None) for every golden test step selector."""
    out = []
    for t in sorted(repo.glob("modules/**/tests/golden.nf.test")) + sorted(repo.glob("subworkflows/**/tests/golden.nf.test")):
        text = t.read_text(encoding="utf-8")
        if 'tag "golden"' not in text:
            raise ValueError(f"{t} is not tagged golden")
        sel = set(SELECT_RE.findall(text))
        if not sel:
            raise ValueError(f"{t} selects no manifest step with it.cmd_name == '...'")
        out += [(str(t.relative_to(repo)), c, w or None) for c, w in sorted(sel)]
    return out


def coverage(manifest: dict, tests: list[tuple[str, str, str | None]]) -> list[dict]:
    rows = []
    for s in manifest["steps"]:
        hits = [t for t, c, w in tests if c == s["cmd_name"] and (w is None or w == s["cwd"])]
        rows.append({**{k: s[k] for k in ("step_index", "name", "cmd_name", "cwd")}, "tests": hits})
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("manifest", type=Path)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    a = ap.parse_args(argv)
    rows = coverage(json.loads(a.manifest.read_text(encoding="utf-8")), find_tests(a.repo))
    print("| # | FragPipe step | Tool (cmd_name) | Work dir | Golden test |")
    print("| - | ------------- | --------------- | -------- | ----------- |")
    for r in rows:
        print(f"| {r['step_index']} | {r['name']} | {r['cmd_name']} | `{r['cwd']}` | {', '.join(f'`{t}`' for t in r['tests']) or 'none'} |")
    n = sum(1 for r in rows if r["tests"])
    print(f"\nCovered: {n} of {len(rows)} steps.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
