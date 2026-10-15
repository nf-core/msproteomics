"""Workflow keys FragPipe 24.0 fills itself when the .workflow omits them (absent-key defaults audit).

Headless FragPipe loads only the user workflow (Fragpipe.java:500, :524) and loadUi sets only the
components the map names (utils/SwingUtils.java:859-869, `if (comp == null) continue;`), so an
absent key keeps its panel's constructor value; the MSFragger tab holds the closed-search defaults
(TabMsfragger.initMore loadDefaults(closed), TabMsfragger.java:515-517, :2088-2096).
Radio buttons: UiRadio.fromString is setSelected(Boolean.valueOf(s)) (utils/swing/UiRadio.java:40-42),
and in a ButtonGroup `false` cannot clear the selected button.
FragPipe paths are relative to FragPipe-GUI/src/main/java/org/nesvilab/fragpipe/ at tag 24.0.
Each test builds a minimal workflow with the key ABSENT and asserts FragPipe's value.
"""

import importlib.util
import json
import subprocess
import sys
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

ASSETS = Path(__file__).resolve().parents[5] / "assets"
FRAGGER_FIXTURES = Path(__file__).resolve().parents[5] / "tests" / "fixtures" / "fragger_params"


@pytest.fixture
def load(tmp_path):
    """Parse a workflow given as {key: value} through parse_workflow_file, as the module does."""

    def make(props):
        path = tmp_path / "test.workflow"
        path.write_text("".join(f"{k}={v}\n" for k, v in props.items()))
        return pfw.parse_workflow_file(str(path))

    return make


# --------------------------------------------------------------------------------------
# Part A: msfragger.* keys other panels read from the MSFragger tab
# --------------------------------------------------------------------------------------


class TestMsfraggerFormValues:
    def test_calibrated_mzml_flag_uses_closed_calibrate_mass(self, load):
        """calibrate_mass absent = closed default 2 (fragger_closed.params:12), so
        isWriteCalMzml() && getMassCalibration() > 0 is true (FragpipeRun.java:1612, :1642, :2090)."""
        p = load({"msfragger.write_calibrated_mzml": "true"})
        out = pfw.generate_json_output(p, ["DDA"])
        assert out["msfragger"]["has_calibrated_mzml"] is True
        assert pfw.opair_spectra_kind(p) == "calibrated"

    def test_ionquant_modmasses_use_closed_var_and_fix_mods(self, load):
        """table.var-mods / table.fix-mods absent = closed variable_mod_01/02 and add_C_cysteine
        (fragger_closed.params:68-69, :126) in getVarModMassSet/getFixedModMassSet
        (TabMsfragger.java:1705-1725, FragpipeRun.java:1852-1855)."""
        p = load({"database.decoy-tag": "rev_"})
        assert pfw._parse_mod_masses(p) == ["15.9949", "42.0106", "57.02146"]

    def test_ionquant_modmasses_keep_closed_fix_mods_when_only_var_mods_given(self, load):
        """table.fix-mods absent alone: add_C_cysteine = 57.02146 stays (fragger_closed.params:126)."""
        p = load({"msfragger.table.var-mods": "79.96633,STY,true,3"})
        assert pfw._parse_mod_masses(p) == ["57.02146", "79.96633"]

    def test_ptmprophet_mods_use_closed_var_mods(self, load):
        """PtmProphetPanel.getMSFraggerMods reads the var-mods table (PtmProphetPanel.java:145-206):
        closed 15.9949 M and 42.0106 [^ (fragger_closed.params:68-69); [^ becomes n."""
        p = load({"database.decoy-tag": "rev_"})
        assert pfw._get_ptmprophet_mods_from_msfragger(p) == "M:15.9949,n:42.0106"

    def test_ptmprophet_offset_merge_uses_closed_restrict_deltamass_to(self, load):
        """restrict_deltamass_to absent = all (fragger_closed.params:44), so getOffsetsAsMods gives the
        20 residues (TabMsfragger.java:2262-2268) and the merge appends them (PtmProphetPanel.java:150-166)."""
        p = load({"msfragger.table.var-mods": "79.96633,STY,true,3", "msfragger.mass_offsets": "79.96633"})
        assert pfw._get_ptmprophet_mods_from_msfragger(p) == "STYACDEFGHIKLMNPQRVW:79.96633"

    def test_peptideprophet_reads_enzyme_name_set_by_dropdown(self, load):
        """search_enzyme_name_1 absent with enzyme-dropdown-1=nonspecific: the dropdown listener fills it
        (TabMsfragger.java:699-753), getEnzymeName() is nonspecific (:2143-2145) and
        CmdPeptideProphet adds --nontt --nonmc (CmdPeptideProphet.java:391-404)."""
        p = load({"msfragger.misc.fragger.enzyme-dropdown-1": "nonspecific"})
        assert pfw.peptideprophet_cmd_opts(p) == "--nontt --nonmc"

    def test_peptideprophet_reads_second_enzyme_cut_set_by_dropdown(self, load):
        """search_enzyme_cut_2 absent with enzyme-dropdown-2=argc: the listener sets cut R
        (TabMsfragger.java:699-753), getEnzymeCut2() is not blank (:2146-2148), so --nontt --nonmc
        (CmdPeptideProphet.java:395)."""
        p = load({"msfragger.misc.fragger.enzyme-dropdown-2": "argc"})
        assert pfw.peptideprophet_cmd_opts(p) == "--nontt --nonmc"

    def test_ptmshepherd_massdiff_to_varmod_uses_closed_default(self, load):
        """mass_diff_to_variable_mod absent = closed 0 (fragger_closed.params:36):
        Integer.toString(getMassDiffToVariableMod()) (FragpipeRun.java:2136, TabMsfragger.java:480-482)."""
        p = load({"ptmshepherd.run-shepherd": "true"})
        assert "msfragger_massdiff_to_varmod = 0" in pfw.generate_ptmshepherd_config(p["ptmshepherd"], p).splitlines()

    def test_ptmshepherd_isotope_error_uses_closed_default(self, load):
        """isotope_error absent = closed 0/1/2 (fragger_closed.params:18), put when not blank
        (FragpipeRun.java:2143-2145)."""
        p = load({"ptmshepherd.run-shepherd": "true", "msfragger.mass_diff_to_variable_mod": "0"})
        assert "isotope_error = 0/1/2" in pfw.generate_ptmshepherd_config(p["ptmshepherd"], p).splitlines()

    def test_workflow_without_msfragger_keys_gets_closed_search_configs(self, tmp_path):
        """No msfragger.* key at all: FragPipe runs MSFragger (UiCheck default true, TabMsfragger.java:570)
        with the closed defaults; all four params files equal the ones the parser writes from the
        msfragger.* keys FragPipe 24.0 itself saved for such a workflow (closed_defaults_saved)."""
        saved = (FRAGGER_FIXTURES / "closed_defaults_saved" / "msfragger_keys.workflow").read_text()
        for name, extra in (("keyless", ""), ("saved", saved)):
            (tmp_path / name).mkdir()
            (tmp_path / name / "w.workflow").write_text("database.decoy-tag=rev_\nworkflow.threads=2\n" + extra)
            common = [sys.executable, str(SCRIPT), *TEMPLATE_ARGS, "--workflow", str(tmp_path / name / "w.workflow"), "--data-types", "DDA"]
            subprocess.run(common + ["--output-json", "--outdir", str(tmp_path / name)], check=True, capture_output=True)
            subprocess.run(common + ["--outdir", str(tmp_path / name / "configs")], check=True, capture_output=True)
        msfragger = json.loads((tmp_path / "keyless" / "tool_configs.json").read_text())["msfragger"]
        assert msfragger["run"] is True
        assert msfragger["calibrate_mass"] == {"DDA": "2", "DIA": "2", "GPF-DIA": "2", "DDA+": "2"}
        for params in ("fragger.params", "fragger_dia.params", "fragger_gpfdia.params", "fragger_dda_plus.params"):
            assert (tmp_path / "keyless" / "configs" / params).read_text() == (tmp_path / "saved" / "configs" / params).read_text(), params

    def test_msbooster_thread_refusal_names_the_pipeline_parameter(self):
        """workflow.threads absent: FragPipe takes availableProcessors()-1 (Fragpipe.java:508-511), which
        the pipeline cannot reproduce; the refusal says how to set it (conf/modules.config --threads)."""
        with pytest.raises(ValueError, match="--fragpipe_threads"):
            pfw.generate_msbooster_config({"run-msbooster": "true"}, {}, True)


# --------------------------------------------------------------------------------------
# Part C: validation and quant panels
# --------------------------------------------------------------------------------------


class TestValidationQuantDefaults:
    def test_percolator_cmd_opts_used_when_run_key_absent(self, load):
        """run-percolator absent: the radio starts selected (PercolatorPanel.java:134, PSMValidation.java:81-82)
        and getCmdOpts() is passed (FragpipeRun.java:1641, CmdPercolator.java:152)."""
        p = load({"percolator.cmd-opts": "--only-psms"})
        assert pfw.generate_tool_args(p, "percolator") == "--only-psms"

    def test_percolator_false_cannot_clear_its_radio(self, load):
        """run-percolator=false, run-peptide-prophet absent: setSelected(false) leaves Percolator selected
        (UiRadio.java:40-42 in a ButtonGroup), so Percolator runs."""
        p = load({"percolator.run-percolator": "false"})
        assert (pfw.get_run_flag(p, "percolator"), pfw.get_run_flag(p, "peptideprophet")) == (True, False)

    def test_percolator_and_peptideprophet_both_false_runs_percolator(self, load):
        """Both radios false: neither call changes the group, Percolator stays selected (PSMValidation.java:81-82)."""
        p = load({"percolator.run-percolator": "false", "peptide-prophet.run-peptide-prophet": "false"})
        assert (pfw.get_run_flag(p, "percolator"), pfw.get_run_flag(p, "peptideprophet")) == (True, False)

    def test_freequant_false_cannot_clear_its_radio(self, load):
        """run-freequant=false, run-ionquant absent, MS1 quant on: FreeQuant stays selected
        (QuantPanelLabelfree.java:234 UiRadio(..., true), :259 IonQuant false; isRunFreeQuant :114-116)."""
        p = load({"quantitation.run-label-free-quant": "true", "freequant.run-freequant": "false"})
        assert (pfw.get_run_flag(p, "freequant"), pfw.get_run_flag(p, "ionquant")) == (True, False)

    @pytest.mark.parametrize(
        "workflow, fillers, expected",
        [
            # FragpipeUiState (FragPipe 24.0's own load + save in the licensed FragPipe 24.0 image,
            # 2026-10-12) saved these run flags for workflows with all four twins true
            ("lfq", 0, {"percolator": False, "peptideprophet": True, "freequant": True, "ionquant": False}),
            ("small", 45, {"percolator": False, "peptideprophet": True, "freequant": False, "ionquant": True}),
            ("small", 20, {"percolator": False, "peptideprophet": True, "freequant": True, "ionquant": False}),
            # `fragpipe --headless --dry-run` (FragPipe 24.0, the flags FRAGPIPEHEADLESS passes, 2026-10-12) ran
            # IonQuant with 96 map keys (table 128) and FreeQuant with 97 (table 256), PeptideProphet in both
            ("headless", 83, {"percolator": False, "peptideprophet": True, "freequant": False, "ionquant": True}),
            ("headless", 84, {"percolator": False, "peptideprophet": True, "freequant": True, "ionquant": False}),
        ],
    )
    def test_both_radio_twins_true_load_like_fragpipe(self, tmp_path, workflow, fillers, expected):
        """Both twins true: the ButtonGroup keeps the radio loadUi sets last, in HashMap order
        (Fragpipe.java:879; SwingUtils.valuesSet :861). LFQ-MBR (413 keys, table 1024) and 20 fillers
        (table 64) end on FreeQuant, 45 fillers (table 128) on IonQuant."""
        twins = ["percolator.run-percolator=true", "peptide-prophet.run-peptide-prophet=true",
                 "freequant.run-freequant=true", "ionquant.run-ionquant=true", "quantitation.run-label-free-quant=true"]
        if workflow == "lfq":
            lines = [line for line in (ASSETS / "LFQ-MBR.workflow").read_text().splitlines()
                     if line.split("=", 1)[0] not in {t.split("=", 1)[0] for t in twins}]
        elif workflow == "headless":
            # run.workflow as headlessWorkflowText writes it: the workflow plus database.db-path
            lines = ["protein-prophet.run-protein-prophet=true", "phi-report.run-report=true", "database.decoy-tag=rev_",
                     "workflow.threads=2"] + [f"filler.key{i}=x" for i in range(1, fillers + 1)] + ["database.db-path=/data/db.fasta"]
        else:
            lines = ["database.decoy-tag=rev_", "workflow.threads=2"] + [f"filler.key{i}=x" for i in range(1, fillers + 1)]
        path = tmp_path / "twins.workflow"
        path.write_text("\n".join(lines + twins) + "\n")
        p = pfw.parse_workflow_file(str(path))
        assert {tool: pfw.get_run_flag(p, tool) for tool in expected} == expected

    def test_freequant_tolerances_default_when_section_present(self, load):
        """freequant.rt-tol / mz-tol absent: toMap() has the spinner values, getFreequantOptsAsText()
        gives " --ptw 0.4 --tol 10" (QuantPanelLabelfree.java:154-161, :238-243)."""
        p = load({"quantitation.run-label-free-quant": "true", "freequant.run-freequant": "true"})
        assert pfw.generate_tool_args(p, "freequant") == "--ptw 0.4 --tol 10"

    def test_freequant_legacy_config_without_tolerances(self, load):
        """The legacy config pass writes freequant.config for any freequant.* section, FreeQuant on or off."""
        p = load({"freequant.run-freequant": "false"})
        assert pfw.generate_tool_config("freequant", p["freequant"], p) == ("freequant.config", "--ptw 0.4 --tol 10")

    def test_freequant_options_without_any_freequant_key(self, load):
        """No freequant.* key, MS1 quant on: FreeQuant is the default radio and gets --ptw 0.4 --tol 10
        (QuantPanelLabelfree.java:154-161, FragpipeRun.java:1817)."""
        p = load({"quantitation.run-label-free-quant": "true"})
        out = pfw.generate_json_output(p, ["DDA"])
        assert (out["freequant"]["run"], out["freequant"]["args"]) == (True, "--ptw 0.4 --tol 10")

    def test_freequant_present_tolerances_take_the_spinner_format(self, load):
        """Present values are written as the spinners format them: "0.0#" and "0.#" (QuantPanelLabelfree.java:238-243)."""
        p = load({"freequant.run-freequant": "true", "freequant.rt-tol": "0.40", "freequant.mz-tol": "10.0"})
        assert pfw.generate_freequant_config(p["freequant"]) == "--ptw 0.4 --tol 10"

    def test_mbg_ionquant_tolerances_default(self, load):
        """ionquant.mztol/rttol/imtol absent with MBG: quantPanelLabelfree.toMap() widget defaults 10, 0.4, 0.05
        (CmdMBGMatch.java:189-195, QuantPanelLabelfree.java:280-285)."""
        mbg = {"max_glycan_q": "0.05", "min_psms": "5", "min_glycans": "3", "fdr": "0.01", "expand_db": "0",
               "max_skips": "1", "allow_chimeric": "false", "residues_to_add": ""}
        p = load({"mbg.run-mbg": "true", **{f"mbg.{k}": v for k, v in mbg.items()}})
        argv = pfw.generate_mbg_args(p["mbg"], p)
        got = {flag: argv[argv.index(flag) + 1] for flag in ("--mztol", "--rttol", "--imtol")}
        assert got == {"--mztol": "10", "--rttol": "0.4", "--imtol": "0.05"}

    def test_fpop_tmt_prot_norm_default_none(self, load):
        """tmtintegrator.prot_norm absent: the combo's first item None (TmtiConfProps.java:78-82,
        TmtiPanel.java:356) gives abundance_*_None.tsv (CmdFpopQuant.java:56-59, :114-145)."""
        p = load({"tmtintegrator.run-tmtintegrator": "true", "fpop.fragpipe.fpop.run-fpop": "true",
                  "fpop.fragpipe.fpop.fpop-tmt": "true"})
        argv = pfw.generate_fpop_config(p["fpop"], p)
        assert (argv[0], argv[-1]) == (
            "tmt-report/abundance_multi-site_None.tsv",
            "tmt-report-unmod/abundance_peptide_None.tsv",
        )

    def test_abacus_label_default_tmt0(self, load):
        """tmtintegrator.channel_num absent: the label combo's first item TMT-0 (TmtiPanel.java:323-324,
        QuantLabel.java:50), which Abacus --plex takes (FragpipeRun.java:1791-1797)."""
        p = load({"tmtintegrator.run-tmtintegrator": "true", "tmtintegrator.extraction_tool": "Philosopher"})
        assert pfw.abacus_iprophet_entries(p, {"filter": {"run": True}})["abacus"]["labels_channel_num"] == "TMT-0"

    def test_tmtintegrator_ms1_int_false_cannot_clear_ms1_radio(self, load):
        """ms1_int=false: the MS1 radio (UiRadio("MS1", null, true)) is the only named button of its group
        (TmtiPanel.java:533-540), so false is a no-op and FragPipe writes ms1_int: true."""
        p = load({"tmtintegrator.run-tmtintegrator": "true", "tmtintegrator.ms1_int": "false"})
        assert "  ms1_int: true" in pfw.generate_tmtintegrator_config(p["tmtintegrator"]).splitlines()
