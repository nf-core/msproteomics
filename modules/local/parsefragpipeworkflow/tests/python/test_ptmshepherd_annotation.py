"""Unit tests for PTM-Shepherd's custom annotation file in parse_fragpipe_workflow.py.

Expected values are FragPipe 24.0 PtmshepherdPanel.toPtmsParamsMap (PtmshepherdPanel.java:216-245):
the annotation_file text field reaches shepherd.config only with the "annotation-glyco" or
"annotation-custom" radio button on; "annotation-unimod" / "annotation-common" write the keywords.
"""

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "resources" / "usr" / "bin" / "parse_fragpipe_workflow.py"
spec = importlib.util.spec_from_file_location("parse_fragpipe_workflow", SCRIPT)
pfw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pfw)

MSFRAGGER = {"mass_diff_to_variable_mod": "0"}


def config_lines(radios, annotation_file):
    params = {"run-shepherd": "true", "annotation_file": annotation_file, **radios}
    text = pfw.generate_ptmshepherd_config(params, {"msfragger": MSFRAGGER, "ptmshepherd": params})
    return dict(line.split(" = ", 1) for line in text.splitlines())


@pytest.mark.parametrize("radio", ["annotation-custom", "annotation-glyco"])
def test_text_field_path_becomes_staging_token(radio):
    cfg = config_lines({radio: "true"}, "/data/my_annotations.txt")
    assert cfg["annotation_file"] == pfw.PTMS_ANNOTATION_FILE_TOKEN == "@ANNOTATION_FILE@"


@pytest.mark.parametrize("radio,keyword", [("annotation-unimod", "unimod"), ("annotation-common", "common")])
def test_keyword_radios_ignore_the_text_field(radio, keyword):
    cfg = config_lines({radio: "true"}, "/data/my_annotations.txt")
    assert cfg["annotation_file"] == keyword


def test_empty_text_field_stays_empty():
    assert config_lines({"annotation-custom": "true"}, "")["annotation_file"] == ""


@pytest.mark.parametrize(
    "radios,text,expected",
    [
        ({"annotation-custom": "true"}, "/data/a.txt", "/data/a.txt"),
        ({"annotation-glyco": "true"}, " /data/a.txt ", "/data/a.txt"),
        ({"annotation-unimod": "true"}, "/data/a.txt", ""),
        ({"annotation-custom": "false", "annotation-common": "true"}, "/data/a.txt", ""),
        ({"annotation-custom": "true"}, "", ""),
        ({"annotation-custom": "true"}, "glyco", ""),
    ],
)
def test_annotation_file_entry(radios, text, expected):
    params = {"annotation_file": text, **radios}
    assert pfw.ptms_annotation_file(params, {"ptmshepherd": params}) == expected
