"""parse_fragpipe_workflow.py MSFragger params file names: FragPipe saves the DDA params as fragger_dda.params when
the manifest also has DIA, GPF-DIA or DIA-Lib files, else fragger.params (CmdMsfragger.java:444-447)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "modules/local/parsefragpipeworkflow/resources/usr/bin/parse_fragpipe_workflow.py"
FIXTURES = REPO / "tests/fixtures/fragpipe_workflows"
sys.path.insert(0, str(SCRIPT.parent))

import parse_fragpipe_workflow as pfw  # noqa: E402

DIA_FILES = {"fragger_dda.params", "fragger_dia.params", "fragger_gpfdia.params", "fragger_dda_plus.params"}
DDA_FILES = {"fragger.params", "fragger_dia.params", "fragger_gpfdia.params", "fragger_dda_plus.params"}


@pytest.mark.parametrize(
    "data_types, name",
    [
        (["DDA"], "fragger.params"),
        (["DDA", "DDA+", "DIA-Quant"], "fragger.params"),
        (["DIA"], "fragger_dda.params"),
        (["DDA", "GPF-DIA"], "fragger_dda.params"),
        (["DIA-Lib"], "fragger_dda.params"),
    ],
)
def test_dda_params_name_follows_the_manifest_data_types(data_types, name):
    assert pfw.msfragger_params_file_name("DDA", data_types) == name


def test_other_data_types_keep_their_names():
    names = {t: pfw.msfragger_params_file_name(t, ["DIA"]) for t in ("DIA", "GPF-DIA", "DDA+")}
    assert names == {"DIA": "fragger_dia.params", "GPF-DIA": "fragger_gpfdia.params", "DDA+": "fragger_dda_plus.params"}


@pytest.mark.parametrize("data_types, files", [("DIA", DIA_FILES), ("DDA", DDA_FILES)])
def test_parser_writes_the_fragpipe_names(tmp_path, data_types, files):
    subprocess.run(
        [sys.executable, str(SCRIPT), "--workflow", str(FIXTURES / "DIA_SpecLib_Quant_diaPASEF.workflow"), "--data-types", data_types,
         "--outdir", str(tmp_path), "--threads", "4"],
        check=True, capture_output=True,
    )
    assert {p.name for p in tmp_path.glob("fragger*.params")} == files


def test_dda_and_fragger_dda_params_have_the_same_content(tmp_path):
    for data_types in ("DIA", "DDA"):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--workflow", str(FIXTURES / "DIA_SpecLib_Quant_diaPASEF.workflow"), "--data-types", data_types,
             "--outdir", str(tmp_path / data_types), "--threads", "4"],
            check=True, capture_output=True,
        )
    assert (tmp_path / "DIA/fragger_dda.params").read_bytes() == (tmp_path / "DDA/fragger.params").read_bytes()
