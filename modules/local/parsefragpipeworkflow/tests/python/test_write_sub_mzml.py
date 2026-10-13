"""Unit tests for Write sub mzML (tab-run.write_sub_mzml) in parse_fragpipe_workflow.py.

Expected values are FragPipe 24.0 (FP = FragPipe-GUI/src/main/java/org/nesvilab/fragpipe/):
TabRun.java:220-243 (checkbox, threshold spinner, Print decoys listener), FragpipeRun.java:1481
(write_mzbin_all), :2388-2403 (CmdWriteSubMzml run and checks), CmdWriteSubMzml.java:72-80.
"""

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "resources" / "usr" / "bin" / "parse_fragpipe_workflow.py"
spec = importlib.util.spec_from_file_location("parse_fragpipe_workflow", SCRIPT)
pfw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pfw)

LFQ_MBR = Path(__file__).resolve().parents[5] / "assets" / "LFQ-MBR.workflow"


def load(tmp_path, replacements=(), drop=()):
    text = LFQ_MBR.read_text()
    for old, new in replacements:
        assert old in text, old
        text = text.replace(old, new)
    text = "\n".join(line for line in text.splitlines() if not any(line.startswith(d) for d in drop))
    wf = tmp_path / "test.workflow"
    wf.write_text(text)
    all_params = pfw.parse_workflow_file(str(wf))
    all_params.setdefault("workflow", {})["threads"] = "4"  # main(): --threads (Fragpipe.java:508-514)
    return all_params


ON = (("tab-run.write_sub_mzml=false", "tab-run.write_sub_mzml=true"),)


def configs(all_params, data_types=("DDA",)):
    pfw.apply_manifest_data_types(all_params, list(data_types))
    pfw.apply_run_tab_listeners(all_params)
    return pfw.generate_json_output(all_params, list(data_types))


class TestRunFlag:
    def test_write_sub_mzml_on_runs_writesubmzml(self, tmp_path):
        out = configs(load(tmp_path, ON))
        assert out["writesubmzml"]["run"] is True
        assert out["writesubmzml"]["prob_threshold"] == "0.5"

    def test_write_sub_mzml_off_does_not_run(self, tmp_path):
        assert configs(load(tmp_path))["writesubmzml"]["run"] is False

    def test_absent_key_keeps_checkbox_default_false(self, tmp_path):
        assert configs(load(tmp_path, drop=("tab-run.write_sub_mzml",)))["writesubmzml"]["run"] is False

    def test_threshold_is_workflow_text(self, tmp_path):
        out = configs(load(tmp_path, ON + (("tab-run.sub_mzml_prob_threshold=0.5", "tab-run.sub_mzml_prob_threshold=0.95"),)))
        assert out["writesubmzml"]["prob_threshold"] == "0.95"

    def test_absent_threshold_is_spinner_initial_value(self, tmp_path):
        out = configs(load(tmp_path, ON, drop=("tab-run.sub_mzml_prob_threshold",)))
        assert out["writesubmzml"]["prob_threshold"] == "0.5"


class TestPrintDecoysListener:
    def test_write_sub_mzml_turns_on_report_decoys(self, tmp_path):
        all_params = load(tmp_path, ON)
        assert all_params["phi-report"]["print-decoys"] == "false"
        out = configs(all_params)
        assert out["filter"]["report_args"].split()[0] == "--decoys"
        assert "--decoys" in out["filter"]["report_args_single_group"].split()

    def test_off_keeps_workflow_print_decoys(self, tmp_path):
        out = configs(load(tmp_path))
        assert "--decoys" not in out["filter"]["report_args"].split()


class TestMzbinAll:
    def test_write_mzbin_all_follows_write_sub_mzml(self, tmp_path):
        on = pfw.generate_msfragger_configs(load(tmp_path, ON)["msfragger"], load(tmp_path, ON))
        off_params = load(tmp_path)
        off = pfw.generate_msfragger_configs(off_params["msfragger"], off_params)
        assert pfw._params_value(on["fragger.params"], "write_mzbin_all") == "1"
        assert pfw._params_value(off["fragger.params"], "write_mzbin_all") == "0"


class TestChecks:
    def test_dda_passes(self, tmp_path):
        pfw.check_write_sub_mzml(load(tmp_path, ON), ["DDA"])

    @pytest.mark.parametrize("data_type", ["DIA", "GPF-DIA", "DIA-Lib", "DDA+"])
    def test_non_dda_data_type_stops(self, tmp_path, data_type):
        with pytest.raises(SystemExit, match="there are non-DDA data types"):
            pfw.check_write_sub_mzml(load(tmp_path, ON), ["DDA", data_type])

    def test_dia_quant_is_not_a_checked_type(self, tmp_path):
        pfw.check_write_sub_mzml(load(tmp_path, ON), ["DDA", "DIA-Quant"])

    def test_msfragger_off_stops(self, tmp_path):
        params = load(tmp_path, ON + (("msfragger.run-msfragger=true", "msfragger.run-msfragger=false"),))
        with pytest.raises(SystemExit, match="Write sub mzML was enabled but MSFragger was not"):
            pfw.check_write_sub_mzml(params, ["DDA"])

    def test_split_database_stops(self, tmp_path):
        params = load(tmp_path, ON + (("msfragger.misc.slice-db=1", "msfragger.misc.slice-db=2"),))
        with pytest.raises(SystemExit, match="split database"):
            pfw.check_write_sub_mzml(params, ["DDA"])

    def test_off_never_stops(self, tmp_path):
        pfw.check_write_sub_mzml(load(tmp_path), ["DIA"])
