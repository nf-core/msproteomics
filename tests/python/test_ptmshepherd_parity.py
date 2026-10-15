"""bin/plan_to_tool_configs.py: PTM-Shepherd's shepherd.config and the PairScans, O-Pair and MBG arguments equal FragPipe 24.0's.

Every expected value is FragPipe's own output: tests/fixtures/ptmshepherd_parity/expected/<case>.json was
captured from a headless `--dry-run --export-plan` run of each case workflow (capture.py, capture.sh).
The cases are the bundled workflows that run PTM-Shepherd or O-Pair plus variants that set a key to a
non-default value, an unformatted number, an invalid combo item, or leave it out (cases.json).
plans/<case>/ is FRAGPIPEPLAN's dry-run work dir of the same case (generate_specs.py, tests/lib/record_fragpipe_plans.nf.test);
a case FragPipe refuses has plans/<case>.refusal.txt instead: the task's exit status and FragPipe's error lines.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
FIXTURES = REPO / "tests/fixtures/ptmshepherd_parity"
sys.path.insert(0, str(FIXTURES))
from capture import cases  # noqa: E402

spec = importlib.util.spec_from_file_location("plan_to_tool_configs", REPO / "bin/plan_to_tool_configs.py")
ptc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ptc)

CASES = cases()


@pytest.mark.parametrize("name", sorted(CASES))
def test_planner_matches_fragpipe(name):
    expected = json.loads((FIXTURES / "expected" / f"{name}.json").read_text())
    if "refusal" in expected:
        refusal = (FIXTURES / "plans" / f"{name}.refusal.txt").read_text()
        assert refusal.splitlines()[0] != "exit status 0" and expected["refusal"] in refusal, refusal
        return
    configs, _files = ptc.tool_configs(ptc.load_workdir(FIXTURES / "plans" / name))
    ptms = configs["ptmshepherd"]["args"].splitlines() if configs["ptmshepherd"]["run"] else None
    argv = {tool: configs[tool]["args"] if configs[tool]["run"] else None for tool in ("pairscans", "mbg", "opair")}
    assert (ptms, argv) == (expected["shepherd_config"], {tool: expected[tool] for tool in argv})


def test_every_case_has_a_fragpipe_capture_and_a_plan():
    assert sorted(p.stem for p in (FIXTURES / "expected").glob("*.json")) == sorted(CASES)
    assert sorted(p.name.removesuffix(".plan_inputs.json") for p in (FIXTURES / "plans").glob("*.plan_inputs.json")) == sorted(CASES)
