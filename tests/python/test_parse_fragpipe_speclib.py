"""parse_fragpipe_workflow.py speclibgen.argv_parts: the values SPECLIBGEN's argv builder (speclibgenCommand)
puts into gen_con_spec_lib.py's argv (CmdSpecLibGen.java:137-202), the same as the SPECLIB_* assignments."""

from __future__ import annotations

import json
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "modules/local/parsefragpipeworkflow/resources/usr/bin/parse_fragpipe_workflow.py"
FIXTURES = REPO / "tests/fixtures/fragpipe_workflows"
KEYS = {
    "rt_cal": "SPECLIB_RT_CAL",
    "im_cal": "SPECLIB_IM_CAL",
    "convert_args": "SPECLIB_CONVERT_ARGS",
    "library_args": "SPECLIB_LIBRARY_ARGS",
    "intermediate": "SPECLIB_INTERMEDIATE",
}


def speclibgen_entry(tmp_path: Path, workflow: str) -> dict:
    subprocess.run(
        [sys.executable, str(SCRIPT), "--workflow", str(FIXTURES / f"{workflow}.workflow"), "--data-types", "DIA",
         "--output-json", "--outdir", str(tmp_path), "--threads", "4"],
        check=True, capture_output=True,
    )
    return json.loads((tmp_path / "tool_configs.json").read_text())["speclibgen"]


@pytest.mark.parametrize("workflow", ["DIA_SpecLib_Quant", "glyco-N-DIA"])
def test_argv_parts_equal_the_sourced_assignments(tmp_path, workflow):
    entry = speclibgen_entry(tmp_path, workflow)
    assigned = dict(line.split("=", 1) for line in entry["args"].splitlines())
    sourced = {key: " ".join(shlex.split(assigned[var])) if assigned[var] not in ("", "''") else "" for key, var in KEYS.items()}
    assert {key: entry["argv_parts"][key] for key in KEYS} == sourced


def test_glyco_psm_conversion_keeps_labile_mode(tmp_path):
    parts = speclibgen_entry(tmp_path, "glyco-N-DIA")["argv_parts"]
    assert "--decoy_prefix rev_" in parts["convert_args"] and "--labile_mods nglyc+" in parts["convert_args"]
