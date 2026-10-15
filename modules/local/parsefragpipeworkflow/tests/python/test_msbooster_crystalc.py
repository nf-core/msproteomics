"""Unit tests for the MSBooster and Crystal-C params of parse_fragpipe_workflow.py.

Expected values are FragPipe 24.0: CmdMSBooster.configure (release jar), MSBoosterPanel model
tables, CrystalcPanel.toParams / crystalc.params, and the msbooster_params.txt FragPipe wrote for
assets/LFQ-MBR.workflow's MSBooster settings (golden raw run, before MSBooster rewrites it).
"""

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "resources" / "usr" / "bin" / "parse_fragpipe_workflow.py"
# FragPipe 24.0 templates the module stages (the script takes them as arguments)
TEMPLATES = Path(__file__).resolve().parents[2] / "assets"
TEMPLATE_ARGS = ["--fragger-template", str(TEMPLATES / "fragger_closed.params"), "--diaumpire-template", str(TEMPLATES / "diaumpire_se.params")]
spec = importlib.util.spec_from_file_location("parse_fragpipe_workflow", SCRIPT)
pfw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pfw)
pfw.FRAGGER_TEMPLATE = TEMPLATES / "fragger_closed.params"
pfw.DIAUMPIRE_TEMPLATE = TEMPLATES / "diaumpire_se.params"

LFQ_MBR_MSBOOSTER = {
    "find-best-im-model": "false",
    "find-best-rt-model": "false",
    "find-best-spectra-model": "false",
    "fragmentation-type": "0",
    "im-model": "DIA-NN",
    "koina-url": "",
    "predict-im": "true",
    "predict-rt": "true",
    "predict-spectra": "true",
    "rt-model": "DIA-NN",
    "run-msbooster": "true",
    "spectra-model": "DIA-NN",
    "spectral-library-path": "",
}

THREADS48 = {"workflow": {"threads": "48"}}


def as_dict(text):
    return dict(line.split(" = ", 1) if " = " in line else (line.rstrip(" ="), "") for line in text.splitlines())


class TestMsboosterConfig:
    def test_lfq_mbr_matches_fragpipe_keys(self):
        cfg = as_dict(pfw.generate_msbooster_config(LFQ_MBR_MSBOOSTER, THREADS48, True))
        assert cfg == {
            "useDetect": "false",
            "renamePin": "1",
            "useRT": "true",
            "useSpectra": "true",
            "useIM": "true",  # the MSBOOSTER module ANDs it with the run's timsTOF flag
            "deletePreds": "false",
            "rtModel": "DIA-NN",
            "spectraModel": "DIA-NN",
            "imModel": "DIA-NN",
            "findBestRtModel": "false",
            "findBestSpectraModel": "false",
            "findBestImModel": "false",
            "KoinaURL": "",
            "rtSearchModelsString": "DIA-NN,Deeplc_hela_hf,AlphaPept_rt_generic,Prosit_2019_irt,Prosit_2020_irt_TMT,Prosit_2024_irt_cit",
            "ms2SearchModelsString": "DIA-NN,ms2pip_2021_hcd,ms2pip_timsTOF2024,ms2pip_ttof5600,ms2pip_immuno_hcd,ms2pip_CID_TMT,ms2pip_iTRAQphospho,AlphaPept_ms2_generic,Prosit_2019_intensity,Prosit_2023_intensity_timsTOF,Prosit_2020_intensity_CID,Prosit_2020_intensity_TMT,Prosit_2020_intensity_HCD,UniSpec,Prosit_2024_intensity_cit,Prosit_2025_intensity_MultiFrag",
            "imSearchModelsString": "DIA-NN,AlphaPept_ccs_generic,IM2Deep",
            "FragmentationType": "auto",
            "numThreads": "48",
        }

    def test_no_corr_keys(self):
        params = dict(LFQ_MBR_MSBOOSTER, **{"use-corr": "true", "use-correlated-features": "true"})
        text = pfw.generate_msbooster_config(params, THREADS48, True)
        assert "useCorr" not in text and "useMultipleCorrelatedFeatures" not in text

    def test_library_path_sets_prediction_files(self):
        cfg = as_dict(pfw.generate_msbooster_config(dict(LFQ_MBR_MSBOOSTER, **{"spectral-library-path": "/lib/a.speclib"}), THREADS48, True))
        assert cfg["spectraPredFile"] == cfg["RTPredFile"] == cfg["IMPredFile"] == "{SPECTRAL_LIBRARY}"

    def test_library_path_is_a_tool_config_field(self):
        # msbooster.spectral-library-path goes to MSBOOSTER as a staged file (entry.spectral_library)
        all_params = {"msbooster": dict(LFQ_MBR_MSBOOSTER, **{"spectral-library-path": "/lib/a.speclib"}), "workflow": {"threads": "48"}}
        assert pfw.generate_json_output(all_params, ["DDA"])["msbooster"]["spectral_library"] == "/lib/a.speclib"
        all_params["msbooster"]["spectral-library-path"] = ""
        assert pfw.generate_json_output(all_params, ["DDA"])["msbooster"]["spectral_library"] == ""

    def test_model_display_names_map_to_keys(self):
        cfg = as_dict(pfw.generate_msbooster_config(dict(LFQ_MBR_MSBOOSTER, **{"rt-model": "DeepLC HeLa HF", "spectra-model": "MS2PIP 2021 HCD"}), THREADS48, True))
        assert (cfg["rtModel"], cfg["spectraModel"]) == ("Deeplc_hela_hf", "ms2pip_2021_hcd")

    def test_model_table_lists_every_named_model(self):
        # ModelsTable.fromString re-enables every row (MSBoosterPanel.f)
        cfg = as_dict(pfw.generate_msbooster_config(dict(LFQ_MBR_MSBOOSTER, **{"table.rt-models": "DIA-NN,true;IM2Deep,false"}), THREADS48, True))
        assert cfg["rtSearchModelsString"] == "DIA-NN,IM2Deep"

    def test_missing_threads_fails_when_msbooster_runs(self):
        with pytest.raises(ValueError, match="needs the thread count"):
            pfw.generate_msbooster_config(LFQ_MBR_MSBOOSTER, {"workflow": {}}, True)

    def test_missing_threads_ok_when_msbooster_off(self):
        cfg = as_dict(pfw.generate_msbooster_config(dict(LFQ_MBR_MSBOOSTER, **{"run-msbooster": "false"}), {"workflow": {}}, False))
        assert "numThreads" not in cfg

    def test_missing_threads_ok_when_validation_tab_off(self):
        # TabValidation.java:58-59, :94: "Run validation tools" off disables MSBooster even with run-msbooster=true
        all_params = {"_global": {"run-validation-tab": "false"}, "workflow": {}, "msbooster": dict(LFQ_MBR_MSBOOSTER)}
        run = pfw.get_run_flag(all_params, "msbooster")
        assert run is False
        cfg = as_dict(pfw.generate_msbooster_config(all_params["msbooster"], all_params, run))
        assert "numThreads" not in cfg

    def test_unknown_model_fails(self):
        with pytest.raises(ValueError, match="not an MSBooster model"):
            pfw.generate_msbooster_config(dict(LFQ_MBR_MSBOOSTER, **{"rt-model": "NoSuchModel"}), THREADS48, True)


class TestCrystalcConfig:
    def test_jar_defaults_regardless_of_workflow(self):
        text = pfw.generate_crystalc_config({"run-crystalc": "true", "precursor_mass": "5"})
        assert text.splitlines() == [
            "precursor_charge = 1 6",
            "isotope_number = 3",
            "precursor_mass = 20",
            "precursor_isolation_window = 0.7",
            "correct_isotope_error = false",
        ]
