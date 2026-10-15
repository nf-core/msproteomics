"""parse_fragpipe_workflow.py: PTM-Shepherd's shepherd.config and the PairScans, O-Pair and MBG arguments equal FragPipe 24.0's.

Every expected value is FragPipe's own output: tests/fixtures/ptmshepherd_parity/expected/<case>.json was
captured from a headless `--dry-run --export-plan` run of each case workflow (capture.py, capture.sh).
The cases are the bundled workflows that run PTM-Shepherd or O-Pair plus variants that set a key to a
non-default value, an unformatted number, an invalid combo item, or leave it out (cases.json).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PARSER = REPO / "modules/local/parsefragpipeworkflow/resources/usr/bin/parse_fragpipe_workflow.py"
# FragPipe 24.0 templates the module stages (the script takes them as arguments)
TEMPLATES = Path(__file__).resolve().parents[2] / "modules" / "local" / "parsefragpipeworkflow" / "assets"
TEMPLATE_ARGS = ["--fragger-template", str(TEMPLATES / "fragger_closed.params"), "--diaumpire-template", str(TEMPLATES / "diaumpire_se.params")]
FIXTURES = REPO / "tests/fixtures/ptmshepherd_parity"
sys.path.insert(0, str(FIXTURES))
from capture import case_workflow, cases  # noqa: E402

CASES = cases()


def parse(name: str, tmp_path: Path) -> subprocess.CompletedProcess:
    workflow = tmp_path / f"{name}.workflow"
    workflow.write_text(case_workflow(CASES[name]))
    return subprocess.run(
        [sys.executable, str(PARSER), *TEMPLATE_ARGS, "--workflow", str(workflow), "--data-types", "DDA", "--output-json",
         "--outdir", str(tmp_path), "--threads", "4"],  # FragPipe ran with --threads 4 (capture.sh)
        capture_output=True, text=True,
    )


@pytest.mark.parametrize("name", sorted(CASES))
def test_parser_matches_fragpipe_plan(name, tmp_path):
    expected = json.loads((FIXTURES / "expected" / f"{name}.json").read_text())
    result = parse(name, tmp_path)
    if "refusal" in expected:
        assert result.returncode != 0 and expected["refusal"] in result.stderr + result.stdout, result.stderr
        return
    assert result.returncode == 0, result.stderr
    configs = json.loads((tmp_path / "tool_configs.json").read_text())
    ptms = configs["ptmshepherd"]["args"].splitlines() if configs["ptmshepherd"]["run"] else None
    argv = {tool: configs[tool]["args"] if configs[tool]["run"] else None for tool in ("pairscans", "mbg", "opair")}
    assert (ptms, argv) == (expected["shepherd_config"], {tool: expected[tool] for tool in argv})


def test_every_case_has_a_fragpipe_capture():
    assert sorted(p.stem for p in (FIXTURES / "expected").glob("*.json")) == sorted(CASES)
