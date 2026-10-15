"""Unit tests for the checks after which FragPipe 24.0 aborts and the saved-form handling of parse_fragpipe_workflow.py.

Each refusal quotes the FragPipe 24.0 source it follows (FragPipe-GUI/src/main/java/org/nesvilab/fragpipe/):
CmdMSBooster.java:106-131, FragpipeRun.java:252-266, CmdMsfragger.java:331-335 and :418-427,
tools/ptmshepherd/PtmshepherdPanel.java:816-829, cmd/CmdIonquant.java:303-304.
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

FIXTURES = Path(__file__).resolve().parents[5] / "tests" / "fixtures" / "fragpipe_workflows"
LFQ_MBR = FIXTURES / "LFQ-MBR.workflow"
ASSETS = Path(__file__).resolve().parents[5] / "assets"


def lfq(**changes):
    """LFQ-MBR (MSBooster, Percolator, IonQuant on) with 'prefix.key' = value overrides."""
    params = pfw.parse_workflow_file(str(LFQ_MBR))
    for dotted, value in changes.items():
        prefix, key = dotted.split("__", 1)
        prefix = prefix.replace("_", "-")
        # MSFragger, O-Pair and PTM-Shepherd keys keep their underscores, the other tools' keys are hyphenated
        key = key if prefix in ("msfragger", "opair", "ptmshepherd") else key.replace("_", "-")
        params.setdefault(prefix, {})[key] = value
    return params


def messages(params):
    return pfw.workflow_refusals(params)


def test_bundled_workflows_are_not_refused():
    for path in sorted(FIXTURES.glob("*.workflow")) + sorted(ASSETS.glob("*.workflow")):
        assert messages(pfw.parse_workflow_file(str(path))) == [], path.name


class TestMsbooster:
    def test_open_search_is_refused(self):
        # CmdMSBooster.java:107-111; TabMsfragger.isOpenSearch: Da units, lower < -3, upper > 3
        params = lfq(msfragger__precursor_mass_units="0", msfragger__precursor_mass_lower="-150", msfragger__precursor_mass_upper="500")
        assert any(m.startswith("MSBooster is incompatible with open search. Please disable MSBooster.") for m in messages(params))

    def test_open_search_with_absent_tolerances_uses_the_panel_spinners(self):
        # TabMsfragger.isOpenSearch reads the spinners (TabMsfragger.java:477); absent keys keep the closed
        # defaults -20/20 (fragger_closed.params:4-5), so Da units alone make it an open search
        params = lfq(msfragger__precursor_mass_units="0")
        del params["msfragger"]["precursor_mass_lower"], params["msfragger"]["precursor_mass_upper"]
        assert any(m.startswith("MSBooster is incompatible with open search. Please disable MSBooster.") for m in messages(params))

    def test_ppm_window_is_no_open_search(self):
        params = lfq(msfragger__precursor_mass_units="1", msfragger__precursor_mass_lower="-150", msfragger__precursor_mass_upper="500")
        assert messages(params) == []

    def test_open_search_without_msbooster_is_allowed(self):
        params = lfq(msfragger__precursor_mass_units="0", msfragger__precursor_mass_lower="-150", msfragger__precursor_mass_upper="500",
                     msbooster__run_msbooster="false")
        assert messages(params) == []

    def test_both_predictions_off_is_refused(self):
        # CmdMSBooster.java:113-116
        params = lfq(msbooster__predict_rt="false", msbooster__predict_spectra="false")
        assert any(m.startswith("At least one of <b>Predict RT</b> or <b>Predict spectra</b> must be enabled if you want to run MSBooster.") for m in messages(params))

    def test_koina_model_without_url_is_refused(self):
        # CmdMSBooster.java:121-124
        params = lfq(msbooster__rt_model="Prosit_2019_irt")
        assert any(m.startswith("Koina URL is required for non DIA-NN models.\nPlease go to <b>Validation</b> tab and adjust the settings.") for m in messages(params))

    def test_find_best_without_url_is_refused(self):
        # CmdMSBooster.java:126-129
        params = lfq(msbooster__find_best_rt_model="true")
        assert any(m.startswith("Koina URL is required for <b>Find best RT/spectral/IM model</b>.\nPlease go to <b>Validation</b> tab") for m in messages(params))

    def test_koina_model_with_url_runs(self):
        params = lfq(msbooster__rt_model="Prosit_2019_irt", msbooster__koina_url="https://koina.example/v2/models/")
        assert messages(params) == []


class TestMsbosterDefaults:
    # MSBoosterPanel.java:142-148: Predict RT, Predict spectra and Predict IM start checked; the find-best boxes unchecked
    def test_absent_predict_keys_are_on_in_the_refusal_and_the_config(self):
        params = lfq()
        for key in ("predict-rt", "predict-spectra", "predict-im"):
            params["msbooster"].pop(key, None)
        assert messages(params) == []
        text = pfw.generate_msbooster_config(params["msbooster"], {"workflow": {"threads": "4"}}, True)
        assert {"useRT = true", "useSpectra = true", "useIM = true"} <= set(text.splitlines())

    def test_explicit_false_for_both_is_refused_and_written_as_false(self):
        params = lfq(msbooster__predict_rt="false", msbooster__predict_spectra="false")
        assert any(m.startswith("At least one of <b>Predict RT</b>") for m in messages(params))
        text = pfw.generate_msbooster_config(params["msbooster"], {"workflow": {"threads": "4"}}, True)
        assert {"useRT = false", "useSpectra = false"} <= set(text.splitlines())

    def test_one_prediction_is_enough(self):
        assert messages(lfq(msbooster__predict_rt="false")) == []


class TestMsbosterModelNames:
    # ui.translations:146-170: the saved workflow holds the config value, MSBoosterPanel.modelMap maps the UI value
    @pytest.mark.parametrize("saved, key", [
        ("Deeplc_hela_hf", "Deeplc_hela_hf"),
        ("Prosit_2024_irt_cit", "Prosit_2024_irt_cit"),
        ("DeepLC HeLa HF", "Deeplc_hela_hf"),
        ("DIA-NN", "DIA-NN"),
    ])
    def test_rt_model(self, saved, key):
        text = pfw.generate_msbooster_config({"rt-model": saved, "workflow-threads": "4"}, {"workflow": {"threads": "4"}}, True)
        assert f"rtModel = {key}" in text.splitlines()

    def test_citrullination_models(self):
        params = {"rt-model": "Prosit_2024_irt_cit", "spectra-model": "Prosit_2024_intensity_cit", "im-model": "DIA-NN"}
        text = pfw.generate_msbooster_config(params, {"workflow": {"threads": "4"}}, True)
        assert "spectraModel = Prosit_2024_intensity_cit" in text.splitlines()

    def test_unknown_model_is_an_error(self):
        with pytest.raises(ValueError, match="is not an MSBooster model"):
            pfw.generate_msbooster_config({"rt-model": "nonsense"}, {"workflow": {"threads": "4"}}, True)


class TestOutputFormat:
    def test_percolator_needs_pin_and_pepxml(self):
        # FragpipeRun.java:252-256
        assert any(m.startswith("Percolator was enabled but MSFragger's output formats did not contain pin or pepXML.") for m in messages(lfq(msfragger__output_format="pepXML")))
        assert any(m.startswith("Percolator was enabled but MSFragger's output formats did not contain pin or pepXML.") for m in messages(lfq(msfragger__output_format="PIN")))

    @pytest.mark.parametrize("fmt", ["pepXML_pin", "tsv_pepXML", "pepXML"])
    def test_peptideprophet_needs_pepxml(self, fmt):
        # FragpipeRun.java:260-264
        params = lfq(percolator__run_percolator="false", peptide_prophet__run_peptide_prophet="true", msfragger__output_format=fmt)
        assert messages(params) == []

    def test_peptideprophet_without_pepxml_is_refused(self):
        params = lfq(percolator__run_percolator="false", peptide_prophet__run_peptide_prophet="true", msfragger__output_format="tsv")
        assert any(m.startswith("PeptideProphet was enabled but MSFragger's output formats did not contain pepXML.") for m in messages(params))


class TestMsfragger:
    def test_labile_remainder_masses_need_localize_delta_mass(self):
        # CmdMsfragger.java:418-427
        params = lfq(msfragger__labile_search_mode="labile", msfragger__remainder_fragment_masses="203.0794",
                     msfragger__localize_delta_mass="false")
        assert any(m.startswith("Fragment remainder masses are specified") for m in messages(params))

    def test_saved_standard_is_the_off_mode(self):
        # ui.translations:23-24: "standard" is saved for the combo's "off"
        params = lfq(msfragger__labile_search_mode="standard", msfragger__remainder_fragment_masses="203.0794",
                     msfragger__localize_delta_mass="false")
        assert messages(params) == []

    def test_labile_with_localize_delta_mass_runs(self):
        params = lfq(msfragger__labile_search_mode="labile", msfragger__remainder_fragment_masses="203.0794",
                     msfragger__localize_delta_mass="true")
        assert messages(params) == []

    def test_remainder_masses_without_labile_mode_run(self):
        params = lfq(msfragger__labile_search_mode="off", msfragger__remainder_fragment_masses="203.0794",
                     msfragger__localize_delta_mass="false")
        assert messages(params) == []

    def test_invalid_var_mod_sites_are_refused(self):
        # CmdMsfragger.java:331-335: pattern [A-Znc\[\]\^\*]+, enabled mods only
        params = lfq(**{"msfragger__table.var-mods": "15.9949,m1,true,3"})
        assert any(m.startswith("Invalid variable modification sites: m1.\nPlease check the tooltip or documentations.") for m in messages(params))

    def test_disabled_mod_with_invalid_sites_runs(self):
        params = lfq(**{"msfragger__table.var-mods": "15.9949,m1,false,3; 42.0106,[^,true,1"})
        assert messages(params) == []

    @pytest.mark.parametrize("value", ["0", "False", ""])
    def test_localize_delta_mass_other_than_true_is_off(self, value):
        # TabMsfragger.java:284: Boolean.parseBoolean
        params = lfq(msfragger__labile_search_mode="labile", msfragger__remainder_fragment_masses="203.0794",
                     msfragger__localize_delta_mass=value)
        assert any(m.startswith("Fragment remainder masses are specified") for m in messages(params))

    def test_msfragger_off_skips_its_checks(self):
        # FragpipeRun.java:1462 configures MSFragger only when it runs; :1558 MSBooster needs MSFragger's pepXML
        params = lfq(**{"msfragger__run-msfragger": "false", "msfragger__table.var-mods": "15.9949,m1,true,3",
                        "msfragger__labile_search_mode": "labile", "msfragger__remainder_fragment_masses": "203.0794",
                        "msfragger__localize_delta_mass": "false", "msbooster__predict-rt": "false",
                        "msbooster__predict-spectra": "false"})
        assert messages(params) == []


class TestPtmshepherd:
    @pytest.mark.parametrize("value", ["Hex:162.0528", "Hex:162.0528, HexNAc:203.0794", "a:-1"])
    def test_valid_varmod_masses(self, value):
        assert messages(lfq(ptmshepherd__varmod_masses=value)) == []

    @pytest.mark.parametrize("value", ["Hex", "Hex:", "Hex:1.2.3", "Hex:162.0528 HexNAc:203.0794"])
    def test_malformed_varmod_masses_are_refused(self, value):
        # PtmshepherdPanel.java:816-829, checked even with PTM-Shepherd off (FragpipeRun.java:2113-2122)
        assert any("There are errors in PTM-Shepherd configuration panel on Report tab." in m for m in messages(lfq(ptmshepherd__varmod_masses=value)))


class TestIonquantArguments:
    def test_whitespace_in_a_one_argument_value_is_refused(self):
        # CmdIonquant.java:303-304 adds each value as one argv element; the pipeline splits on spaces
        params = lfq(ionquant__excludemods="M15.9949; STY79.96633")
        assert any(m.startswith("ionquant.excludemods=") for m in messages(params))

    @pytest.mark.parametrize("key", ["light", "medium", "heavy"])
    def test_leftover_label_without_use_labeling_is_skipped(self, key):
        # CmdIonquant.java:296-301: the label strings are skipped unless the use-labeling radio is on
        params = lfq(**{f"ionquant__{key}": "K4.025107; R6.020129", "ionquant__use_labeling": "false"})
        assert messages(params) == []
        params = lfq(**{f"ionquant__{key}": "K4.025107; R6.020129", "ionquant__use_labeling": "true"})
        assert any(m.startswith(f"ionquant.{key}=") for m in messages(params))

    def test_value_without_whitespace_runs(self):
        assert messages(lfq(ionquant__excludemods="M15.9949;STY79.96633")) == []


class TestSpeclibgenRtCal:
    def settings(self, rt_cal, **extra):
        params = {"easypqp.rt-cal": rt_cal, **extra}
        return pfw.speclibgen_settings(params, {})

    def test_saved_tsv_file_form(self):
        # ui.translations:79: "User provided RT calibration file" is saved as "a tsv file"
        settings = self.settings("a tsv file", **{"easypqp.select-file.text": "/x/rt.tsv"})
        assert settings["rt_cal_file"] == "/x/rt.tsv" and settings["rt_cal"] == ""

    def test_ui_text_form(self):
        settings = self.settings("User provided RT calibration file", **{"easypqp.select-file.text": "/x/rt.tsv"})
        assert settings["rt_cal_file"] == "/x/rt.tsv"

    def test_tsv_file_without_a_file_is_an_error(self):
        with pytest.raises(ValueError, match="select-file.text is empty"):
            self.settings("a tsv file")


class TestDiaumpire:
    def test_keys_without_a_panel_component_keep_the_template_value(self):
        # UmpirePanel.collect() overlays only names with a component (UmpirePanel.java:224-244)
        template = pfw.generate_diaumpire_config({})
        written = pfw.generate_diaumpire_config({"WindowType": "TIMSTOF", "WindowSize": "99", "SE.MinMSIntensity": "5", "SE.MinMSMSIntensity": "5"})
        assert written == template

    def test_a_component_key_is_applied(self):
        assert "RPmax = 33" in pfw.generate_diaumpire_config({"RPmax": "33"}).splitlines()

    def test_text_is_the_template_byte_for_byte_with_the_thread_placeholder(self):
        # writePropertiesContent println()s every template line, the final empty line included
        # (PropertiesUtils.java:291-321): FragPipe's umpire-se.params ends "==window setting end\n\n"
        template = pfw.DIAUMPIRE_TEMPLATE.read_text(encoding="utf-8")
        written = pfw.generate_diaumpire_config({})
        assert written.replace("Thread = THREADS", "Thread = 7") == template
        assert written.endswith("==window setting end\n\n")


class TestDiannCmdOpts:
    def test_cmd_opts_is_one_field_not_a_word_of_the_args(self):
        # CmdDiann.java:386-388: cmd.add(additionalCmdOpts), one argv element
        params = {"diann": {"cmd-opts": "--mod 280.1899,280.1899", "run-dia-nn": "true"}}
        assert pfw.diann_settings(params["diann"])["cmd_opts"] == "--mod 280.1899,280.1899"
        assert "--mod" not in pfw.generate_diann_post_args(params["diann"])
