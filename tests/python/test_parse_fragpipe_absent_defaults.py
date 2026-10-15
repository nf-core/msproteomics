"""parse_fragpipe_workflow.py: workflow keys that are absent take FragPipe 24.0's panel defaults.

A headless FragPipe run loads only the user workflow, and SwingUtils.valuesSet sets a component only
when its key is in the file (utils/SwingUtils.java:859-868); an absent key keeps the value the panel
constructor gave the widget. FP = FragPipe-GUI/src/main/java/org/nesvilab/fragpipe/ at tag 24.0.
Every expected value below follows the cited FragPipe line, not the parser's output.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "modules/local/parsefragpipeworkflow/resources/usr/bin"))
import parse_fragpipe_workflow as pfw  # noqa: E402

# FragPipe 24.0 templates the module stages (the script takes them as arguments)
TEMPLATES = REPO / "modules" / "local" / "parsefragpipeworkflow" / "assets"
pfw.FRAGGER_TEMPLATE = TEMPLATES / "fragger_closed.params"
pfw.DIAUMPIRE_TEMPLATE = TEMPLATES / "diaumpire_se.params"

MSFRAGGER = {"msfragger": {"mass_diff_to_variable_mod": "0"}}


def shepherd_lines(params: dict) -> dict:
    """shepherd.config property lines as {key: value}."""
    text = pfw.generate_ptmshepherd_config(params, MSFRAGGER)
    return dict(line.split(" = ", 1) for line in text.splitlines())


class TestPtmShepherdPanelDefaults:

    def test_absent_spectra_cond_ratio_is_panel_default(self):
        # FP/tools/ptmshepherd/PtmshepherdPanel.java:742 UiSpinnerDouble.builder(0.0001, ...), format "0.#######"
        assert shepherd_lines({"run-shepherd": "true"})["spectra_condRatio"] == "0.0001"

    def test_present_spectra_cond_ratio_is_kept(self):
        assert shepherd_lines({"run-shepherd": "true", "spectra_condRatio": "0.01"})["spectra_condRatio"] == "0.01"

    def test_absent_glycan_panel_keys_are_panel_defaults(self):
        # FP/tools/ptmshepherd/PTMSGlycanAssignPanel.java:162 (50.0, "0.#"), :168, :171 (UiSpinnerInt -1 / 3),
        # :199 (0.5, "0.0#")
        lines = shepherd_lines({"run-shepherd": "true", "run_glyco_mode": "true"})
        assert (lines["glyco_ppm_tol"], lines["prob_mass"], lines["glyco_isotope_min"], lines["glyco_isotope_max"]) == (
            "50", "0.5", "-1", "3")

    def test_one_isotope_bound_keeps_the_other_panel_default(self):
        # with only glyco_isotope_max set, FragPipe writes min = -1 too; PTM-Shepherd needs both
        # (GlycoParams.parseGlycoIsotopesParam) or it falls back to -1..3 and ignores the present one
        lines = shepherd_lines({"run-shepherd": "true", "run_glyco_mode": "true", "glyco_isotope_max": "1"})
        assert (lines["glyco_isotope_min"], lines["glyco_isotope_max"]) == ("-1", "1")

    def test_glycan_keys_stay_out_without_glyco_mode(self):
        # getGlycanAssignParams() reaches shepherd.config only when glycan assignment runs (FragpipeRun.java:2137-2138)
        assert "glyco_ppm_tol" not in shepherd_lines({"run-shepherd": "true"})


# FP/tools/opair/OPairPanel.java:118, 122, 127, 132, 135, 150, 153 (and CmdOPair.java:177-209 builds the argv)
OPAIR_DEFAULT_ARGV = [
    "-b 20.0",
    "-c 20.0",
    "-f @TOOLS_DIR@/opair/Glycan_Mods/oxonium_filter.txt",
    "-m 0.05",
    "-x @TOOLS_DIR@/Glycan_Databases/glycan_residues.txt",
    "-y @TOOLS_DIR@/Glycan_Databases/glycan_mods.txt",
    "-n 2",
    "-z ",
    "-t @THREADS@",
    "-i 0",
    "-j 2",
    "-r @FILELIST@",
    "-s @PSM@",
    "-o @OUTDIR@",
]


class TestOPairPanelDefaults:

    def test_absent_keys_give_panel_default_argv(self):
        assert pfw.generate_opair_args({}) == OPAIR_DEFAULT_ARGV

    def test_absent_pairscans_keys_are_panel_defaults(self):
        # OPairPanel.java:139-147: activation combos default to their first item, HCD
        # (tools/enums/ActivationTypes.java:21, UiUtils.createUiCombo(Enum[]) -> Enum::name); both checks false
        assert pfw.generate_pairscans_args({}) == ["HCD", "HCD", "false", "false"]

    def test_hard_stops_pass_with_absent_activation_keys(self):
        # HCD with single_scan_type=false passes CmdOPair.java:89-96, so FragPipe runs
        pfw.check_glyco_hard_stops({"opair": {"run-opair": "true"}})

    def test_pairscans_runs_when_single_scan_type_is_absent(self):
        # FragpipeRun.java:2067-2081 runs PairScans when !isSingleScanType(), and the check defaults to false
        assert pfw._opair_value({}, "single_scan_type") == "false"


class TestMbgPanelDefaults:

    IONQUANT = {"ionquant": {"mztol": "10", "rttol": "0.4", "imtol": "0.05"}}

    def test_absent_keys_give_panel_default_argv(self):
        # FP/tools/mbg/MBGPanel.java:98, 104, 109, 113, 117, 121, 125
        argv = pfw.generate_mbg_args({"residues_to_add": "HexNAc"}, self.IONQUANT)
        flags = dict(zip(argv[1::2], argv[2::2])) if argv[0] == "--match" else {}
        assert {k: flags[k] for k in ("--maxq", "--minpsms", "--minglycans", "--fdr", "--expanddb", "--maxskips", "--allowchimeric")} == {
            "--maxq": "0.01", "--minpsms": "5", "--minglycans": "2", "--fdr": "0.01",
            "--expanddb": "2", "--maxskips": "0", "--allowchimeric": "false"}


class TestPtmProphetOverrideDefault:

    def test_override_without_cmdline_uses_the_text_field_default(self):
        # FP/tools/ptmprophet/PtmProphetPanel.java:86 uiTextCmd = ...text(getDefaults()), getDefaults() reads
        # ptmprophet.cmdline.default (FP/Bundle.properties:46); getCmdLineOpts() returns that text (:241-247)
        assert pfw.generate_ptmprophet_config({"override-defaults": "true"}, {}) == (
            "NOSTACK KEEPOLD STATIC FRAGPPMTOL=10 EM=1 NIONS=b STY:79.966331,M:15.9949 MINPROB=0.5")

    def test_override_with_cmdline_keeps_it(self):
        assert pfw.generate_ptmprophet_config({"override-defaults": "true", "cmdline": "MINPROB=0.9"}, {}) == "MINPROB=0.9"


class TestTransferLearningIsRefused:

    def test_diann_entry_has_no_predicted_library(self):
        # Transfer learning is refused before any output (check_supported_stages, run from main()), so the
        # predicted-library path of CmdDiann.java:131-135, :321-322 and its TransferLearningPanel defaults
        # (predict-spectral-library true, output-format speclib) are not modelled
        all_params = {"diann": {"run-dia-nn": "true"}, "transfer-learning": {"predict-spectral-library": "true"}}
        assert "predicted_library" not in pfw.generate_json_output(all_params, ["DIA"])["diann"]

    def test_transfer_learning_stops_the_parser(self):
        import pytest
        with pytest.raises(SystemExit, match="Transfer learning was enabled"):
            pfw.check_supported_stages({"transfer-learning": {"run-transfer-learning": "true"}})
