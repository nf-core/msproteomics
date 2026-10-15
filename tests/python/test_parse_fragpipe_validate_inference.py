"""parse_fragpipe_workflow.py: validation and inference decisions vs FragPipe 24.0.

FP = FragPipe-GUI/src/main/java/org/nesvilab/fragpipe/ at tag 24.0. Every expected value
below follows the cited FragPipe source line, not the parser's current output.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "modules/local/parsefragpipeworkflow/resources/usr/bin"))
import parse_fragpipe_workflow as pfw  # noqa: E402

# FragPipe 24.0 templates the module stages (the script takes them as arguments)
TEMPLATES = REPO / "modules" / "local" / "parsefragpipeworkflow" / "assets"
pfw.FRAGGER_TEMPLATE = TEMPLATES / "fragger_closed.params"
pfw.DIAUMPIRE_TEMPLATE = TEMPLATES / "diaumpire_se.params"


def parse_with_threads(path: str) -> dict:
    """parse_workflow_file plus the run's thread count: the pipeline passes --threads (workflow.threads), which
    the parser requires whenever MSBooster runs (its RT loess output depends on it)."""
    params = pfw.parse_workflow_file(path)
    params.setdefault("workflow", {})["threads"] = "4"
    return params


FIXTURES = REPO / "tests/fixtures/fragpipe_workflows"


def params(**sections):
    """all_params from {section: {key: value}} (keys without the section prefix)."""
    return {k.replace("_", "-"): v for k, v in sections.items()}


class TestPercolatorMinProb:
    # FP/tools/percolator/PercolatorPanel.java:77-79,140 and FP/cmd/CmdPercolator.java:298
    def test_min_prob_is_java_double_string(self):
        assert pfw.percolator_min_prob(params(percolator={"min-prob": "0.7"})) == "0.7"

    def test_min_prob_default_is_panel_default(self):
        assert pfw.percolator_min_prob({}) == "0.5"

    def test_min_prob_rounds_to_spinner_format(self):
        assert pfw.percolator_min_prob(params(percolator={"min-prob": "0.555"})) == "0.56"

    def test_min_prob_one_prints_as_double(self):
        assert pfw.percolator_min_prob(params(percolator={"min-prob": "1"})) == "1.0"

    @pytest.mark.parametrize("value,text", [(0.5, "0.5"), (1e-4, "1.0E-4"), (1.5e-5, "1.5E-5"), (1.2345678e7, "1.2345678E7")])
    def test_java_double_to_string(self, value, text):
        assert pfw.java_double_to_string(value) == text

    def test_json_carries_min_prob_and_keep_tsv(self):
        # DIA_SpecLib_Quant.workflow:235 percolator.min-prob=0.7
        out = pfw.generate_json_output(parse_with_threads(str(FIXTURES / "DIA_SpecLib_Quant.workflow")), ["DIA"])
        assert (out["percolator"]["min_prob"], out["percolator"]["keep_tsv_files"]) == ("0.7", False)


def edited_workflow(tmp_path, name: str, drop=(), set_=None) -> dict:
    """A fixture workflow with keys removed (drop) or set (set_), parsed with threads."""
    lines = (FIXTURES / name).read_text().splitlines()
    lines = [ln for ln in lines if ln.split("=", 1)[0] not in drop]
    for key, value in (set_ or {}).items():
        lines = [f"{key}={value}" if ln.split("=", 1)[0] == key else ln for ln in lines]
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return parse_with_threads(str(path))


class TestPercolatorRunDefault:
    # FP/tools/percolator/PercolatorPanel.java:134 new UiRadio("Run Percolator", null, true): an absent
    # percolator.run-percolator leaves the radio selected, and CmdPercolator.java:118 (setCmdLineParams) /
    # :152 (addFreeCommandLineParams) pass percolator.cmd-opts. A FragPipe 24.0 --dry-run of
    # LFQ-MBR-paper without the key runs `percolator --only-psms --no-terminate --post-processing-tdc
    # --num-threads ...`, the same argv as with run-percolator=true.
    FRAGPIPE_OPTS = "--only-psms --no-terminate --post-processing-tdc"

    def test_absent_run_key_keeps_cmd_opts(self, tmp_path):
        out = pfw.generate_json_output(edited_workflow(tmp_path, "LFQ-MBR.workflow", drop=("percolator.run-percolator",)), ["DDA"])
        assert (out["percolator"]["run"], out["percolator"]["args"]) == (True, self.FRAGPIPE_OPTS)

    def test_absent_run_key_args_equal_explicit_run(self, tmp_path):
        explicit = pfw.generate_json_output(parse_with_threads(str(FIXTURES / "LFQ-MBR.workflow")), ["DDA"])
        absent = pfw.generate_json_output(edited_workflow(tmp_path, "LFQ-MBR.workflow", drop=("percolator.run-percolator",)), ["DDA"])
        assert absent["percolator"] == explicit["percolator"]

    def test_peptideprophet_selected_drops_percolator_args(self, tmp_path):
        # PSMValidation.java:80-82 one ButtonGroup: PeptideProphet selected deselects Percolator
        all_params = edited_workflow(
            tmp_path, "LFQ-MBR.workflow", drop=("percolator.run-percolator",),
            set_={"peptide-prophet.run-peptide-prophet": "true"},
        )
        out = pfw.generate_json_output(all_params, ["DDA"])
        assert (out["percolator"]["run"], out["percolator"]["args"]) == (False, "")


class TestValidationTab:
    # FP/tabs/TabValidation.java:58-59 "Run validation tools" (run-validation-tab, default true), :94
    # setEnablementUpdater disables Crystal-C, MSBooster, PSM validation, PTMProphet, ProteinProphet and
    # Report; each isRun() is isEnabledAndChecked. A FragPipe 24.0 --dry-run of LFQ-MBR-paper with
    # run-validation-tab=false plans CheckCentroid, MSFragger (+ moves) and IonQuant only.
    GATED = ("crystalc", "msbooster", "percolator", "peptideprophet", "ptmprophet", "proteinprophet", "filter")

    def test_off_disables_validation_tools(self, tmp_path):
        all_params = edited_workflow(
            tmp_path, "LFQ-MBR.workflow",
            set_={"run-validation-tab": "false", "quantitation.run-label-free-quant": "false"},
        )
        assert {t: pfw.get_run_flag(all_params, t) for t in self.GATED + ("msfragger",)} == {
            **{t: False for t in self.GATED}, "msfragger": True}

    def test_on_keeps_validation_tools(self):
        all_params = parse_with_threads(str(FIXTURES / "LFQ-MBR.workflow"))
        assert [t for t in ("msbooster", "percolator", "proteinprophet", "filter") if pfw.get_run_flag(all_params, t)] == [
            "msbooster", "percolator", "proteinprophet", "filter"]

    def test_off_without_psm_readers_is_supported(self, tmp_path):
        all_params = edited_workflow(
            tmp_path, "LFQ-MBR.workflow",
            set_={"run-validation-tab": "false", "quantitation.run-label-free-quant": "false"},
        )
        pfw.check_supported_stages(all_params)

    def test_off_with_ionquant_fails(self, tmp_path):
        # FragPipe plans IonQuant on <group>/psm.tsv (CmdIonquant.java:318), which Report never writes
        all_params = edited_workflow(tmp_path, "LFQ-MBR.workflow", set_={"run-validation-tab": "false"})
        with pytest.raises(SystemExit, match="run-validation-tab=false.*these tools read it: ionquant"):
            pfw.check_supported_stages(all_params)


class TestDigestOnly:
    # FP/tabs/TabMsfragger.java:1932-1934 isRun() = isEnabledAndChecked(checkRun) && !isDigestOnly():
    # a FragPipe 24.0 --dry-run of LFQ-MBR-paper with msfragger.digest-only=true plans no MSFragger step.
    def test_digest_only_fails(self, tmp_path):
        all_params = edited_workflow(tmp_path, "LFQ-MBR.workflow", set_={"msfragger.digest-only": "true"})
        with pytest.raises(SystemExit, match=r"msfragger\.digest-only=true"):
            pfw.check_supported_stages(all_params)

    def test_search_is_supported(self):
        pfw.check_supported_stages(parse_with_threads(str(FIXTURES / "LFQ-MBR.workflow")))


class TestFreequantWithoutReport:
    # FP/FragpipeRun.java:1690 DbAnnotate and :1701 Filter run on reportPanel.isRun(), :1767 Report on
    # reportPanel.isRun() || isRunFreeQuant(). A FragPipe 24.0 --dry-run with FreeQuant and run-report=false
    # plans FreeQuant and PhilosopherReport but no DbAnnotate or Filter.
    FREEQUANT = {"quantitation.run-label-free-quant": "true", "freequant.run-freequant": "true", "ionquant.run-ionquant": "false"}

    def test_report_off_fails(self, tmp_path):
        all_params = edited_workflow(tmp_path, "LFQ-MBR.workflow", set_={**self.FREEQUANT, "phi-report.run-report": "false"})
        with pytest.raises(SystemExit, match="FreeQuant runs without the Philosopher report panel"):
            pfw.check_supported_stages(all_params)

    def test_validation_tab_off_fails(self, tmp_path):
        all_params = edited_workflow(tmp_path, "LFQ-MBR.workflow", set_={**self.FREEQUANT, "run-validation-tab": "false"})
        with pytest.raises(SystemExit, match="FreeQuant runs without the Philosopher report panel"):
            pfw.check_supported_stages(all_params)

    def test_report_off_with_reader_names_only_the_freequant_error(self, tmp_path):
        # FreeQuant runs Report (FragpipeRun.java:1767), so the no-psm.tsv message would be false here
        all_params = edited_workflow(tmp_path, "LFQ-MBR.workflow", set_={**self.FREEQUANT, "phi-report.run-report": "false", "ptmshepherd.run-shepherd": "true"})
        with pytest.raises(SystemExit) as err:
            pfw.check_supported_stages(all_params)
        assert "FreeQuant runs without" in str(err.value) and "no psm.tsv is written" not in str(err.value)

    def test_report_on_is_supported(self, tmp_path):
        pfw.check_supported_stages(edited_workflow(tmp_path, "LFQ-MBR.workflow", set_=self.FREEQUANT))

    def test_ionquant_with_report_off_fails(self, tmp_path):
        # DbAnnotate, Filter and Report all off (FragpipeRun.java:1690, :1701, :1767), but IonQuant reads
        # <group>/psm.tsv (CmdIonquant.java:318)
        all_params = edited_workflow(tmp_path, "LFQ-MBR.workflow", set_={"phi-report.run-report": "false"})
        with pytest.raises(SystemExit, match=r"run-report=false.*these tools read it: ionquant"):
            pfw.check_supported_stages(all_params)


class TestPeptideProphetOptions:
    # FP/cmd/CmdPeptideProphet.java:375-404
    def _opts(self, cmd_opts, enzyme, cut2=""):
        return pfw.peptideprophet_cmd_opts(
            params(peptide_prophet={"cmd-opts": cmd_opts}, msfragger={"search_enzyme_name_1": enzyme, "search_enzyme_cut_2": cut2})
        )

    def test_trypsin_keeps_cmd_opts(self):
        assert self._opts("--decoyprobs --ppm", "stricttrypsin") == "--decoyprobs --ppm"

    def test_nonspecific_adds_nontt_nonmc(self):
        assert self._opts("--ppm", "nonspecific") == "--ppm --nontt --nonmc"

    def test_second_enzyme_cut_adds_nontt_nonmc(self):
        assert self._opts("--ppm", "stricttrypsin", "@") == "--ppm --nontt --nonmc"

    def test_custom_adds_enzyme_nonspecific(self):
        assert self._opts("--ppm", "custom") == "--ppm --nontt --nonmc --enzyme nonspecific"

    def test_nocleavage_keeps_existing_enzyme(self):
        assert self._opts("--enzyme trypsin --nontt", "nocleavage") == "--enzyme trypsin --nontt --nonmc"

    def test_null_enzyme_name_is_not_nocleavage(self):
        assert self._opts("--ppm", "null") == "--ppm"

    def test_semi_specific_named_enzyme_gets_no_flags(self):
        p = params(peptide_prophet={"cmd-opts": "--ppm"}, msfragger={"search_enzyme_name_1": "stricttrypsin", "num_enzyme_termini": "0"})
        assert pfw.peptideprophet_cmd_opts(p) == "--ppm"

    def test_literal_peptideprophet_token_dropped(self):
        assert self._opts("peptideprophet --ppm", "stricttrypsin") == "--ppm"

    # FP/cmd/CmdPeptideProphet.java:190-202
    def test_combine_from_checkbox(self):
        assert pfw.peptideprophet_combine(params(peptide_prophet={"combine-pepxml": "true"}))

    def test_combine_from_cmd_opts(self):
        assert pfw.peptideprophet_combine(params(peptide_prophet={"combine-pepxml": "false", "cmd-opts": "--COMBINE"}))

    def test_open_workflow_combines(self):
        # Open.workflow:231-232
        out = pfw.generate_json_output(parse_with_threads(str(FIXTURES / "Open.workflow")), ["DDA"])
        assert out["peptideprophet"]["combine"] is True


class TestPtmProphetOptions:
    # FP/Bundle.properties:46, FP/tools/ptmprophet/PtmProphetPanel.java:113-138,241-247
    def test_auto_keeps_fragppmtol_10(self):
        p = params(msfragger={"table.var-mods": "15.9949,M,true,3", "fragment_mass_units": "1", "fragment_mass_tolerance": "20"})
        assert pfw.generate_ptmprophet_config({}, p) == "NOSTACK KEEPOLD STATIC FRAGPPMTOL=10 EM=1 NIONS=b M:15.9949 MINPROB=0.5"

    def test_auto_keeps_fragppmtol_for_dalton_units(self):
        p = params(msfragger={"table.var-mods": "15.9949,M,true,3", "fragment_mass_units": "0", "fragment_mass_tolerance": "0.02"})
        assert "FRAGPPMTOL=10" in pfw.generate_ptmprophet_config({}, p).split()

    def test_override_is_verbatim(self):
        assert pfw.generate_ptmprophet_config({"override-defaults": "true", "cmdline": "KEEPOLD STY:79.966 MINPROB=0.9"}, {}) == "KEEPOLD STY:79.966 MINPROB=0.9"

    def test_override_empty_stays_empty(self):
        assert pfw.generate_ptmprophet_config({"override-defaults": "true", "cmdline": ""}, {}) == ""

    # FP/cmd/CmdPtmProphet.java:115-117,138-151
    def test_philosopher_style_is_translated(self):
        cmd = "--keepold --static --em 1 --mods STY:79.966331 --minprob 0.5"
        assert pfw.generate_ptmprophet_config({"override-defaults": "true", "cmdline": cmd}, {}) == "KEEPOLD STATIC EM=1 STY:79.966331 MINPROB=0.5"


class TestFilterOptions:
    # FP/cmd/CmdPhilosopherFilter.java:82-97
    def test_razor_removed(self):
        assert pfw.fragpipe_filter_options("--sequential --razor --prot 0.05", False, False) == "--sequential --prot 0.05"

    def test_ionquant_forces_prot(self):
        assert pfw.fragpipe_filter_options("--sequential --prot 0.05", True, False) == "--sequential --prot 0.01"

    def test_ionquant_appends_prot(self):
        assert pfw.fragpipe_filter_options("--picked", True, False) == "--picked --prot 0.01"

    def test_no_protxml_strips_after_ionquant(self):
        assert pfw.fragpipe_filter_options("--sequential --picked", True, True) == "--picked"

    def test_empty_text(self):
        assert pfw.fragpipe_filter_options("  ", True, False) == ""

    # FP/FragpipeRun.java:1716-1719,1752-1755
    def test_no_proteinprophet_means_no_protxml(self):
        p = params(phi_report={"run-report": "true", "filter": "--sequential --prot 0.01"}, protein_prophet={"run-protein-prophet": "false"})
        entry = pfw.philosopher_filter_entry(p)
        assert (entry["dont_use_protxml"], entry["args"]) == (True, "")

    def test_checkbox_honoured_with_proteinprophet(self):
        p = params(phi_report={"filter": "--sequential --picked --prot 0.01", "dont-use-prot-proph-file": "true"}, protein_prophet={"run-protein-prophet": "true"})
        assert pfw.philosopher_filter_entry(p)["args"] == "--picked"

    # FP/cmd/CmdPhilosopherReport.java:54-62, FP/FragpipeRun.java:1768-1780
    def test_report_msstats_single_group_only(self):
        p = params(
            phi_report={"print-decoys": "true", "remove-contaminants": "true"},
            protein_prophet={"run-protein-prophet": "true"},
            tmtintegrator={"run-tmtintegrator": "true", "philosopher-msstats": "true", "extraction_tool": "Philosopher"},
        )
        entry = pfw.philosopher_filter_entry(p)
        assert (entry["report_args"], entry["report_args_single_group"]) == ("--decoys --removecontam", "--decoys --msstats --removecontam")

    def test_msstats_with_ionquant_extraction_fails(self):
        p = params(tmtintegrator={"run-tmtintegrator": "true", "philosopher-msstats": "true", "extraction_tool": "IonQuant"})
        with pytest.raises(ValueError, match="FragpipeRun.java:1770-1773"):
            pfw.philosopher_filter_entry(p)


class TestAbacusIprophet:
    # FP/FragpipeRun.java:1782-1809, FP/cmd/CmdPhilosopherAbacus.java:106-130, FP/cmd/CmdIprophet.java:44-52
    def _entries(self, name):
        return pfw.generate_json_output(parse_with_threads(str(FIXTURES / f"{name}.workflow")), ["DDA"])

    def test_basic_search_runs_abacus(self):
        # Basic-Search.workflow:239,241,312: prot-level summary, no LFQ
        out = self._entries("Basic-Search")
        assert (out["abacus"]["run"], out["abacus"]["args"], out["iprophet"]["run"]) == (True, "--picked --razor --reprint --tag rev_ --protein", False)

    def test_lfq_ionquant_disables_abacus(self):
        assert self._entries("LFQ-MBR")["abacus"]["run"] is False

    def test_tmt_ionquant_extraction_disables_abacus(self):
        assert self._entries("TMT10")["abacus"]["run"] is False

    def test_peptide_summary_runs_iprophet(self):
        p = params(
            phi_report={"run-report": "true", "pep-level-summary": "true", "prot-level-summary": "false", "filter": "--uniqueonly --picked --uniqueonly"},
            database={"decoy-tag": "XXX_"},
        )
        entries = pfw.abacus_iprophet_entries(p, {"filter": {"run": True}})
        assert (entries["abacus"]["args"], entries["iprophet"]["args"]) == (
            "--uniqueonly --picked --razor --reprint --tag XXX_ --peptide",
            "--decoy XXX_ --nonsp --output combined",
        )

    def test_philosopher_tmt_extraction_passes_label(self):
        p = params(
            phi_report={"run-report": "true"},
            tmtintegrator={"run-tmtintegrator": "true", "extraction_tool": "Philosopher", "channel_num": "TMT-16"},
        )
        assert pfw.abacus_iprophet_entries(p, {"filter": {"run": True}})["abacus"]["labels_channel_num"] == "TMT-16"


class TestSaintexpress:
    # FP/cmd/CmdSaintExpress.java:85-90, FP/tools/downstream/SaintexpressPanel.java:111-115
    def test_args_and_single_cmd_opts(self):
        p = params(saintexpress={"run-saint-express": "true", "max-replicates": "3", "virtual-controls": "50", "cmd-opts": " -f -x "})
        out = pfw.generate_json_output(p, ["DDA"])
        assert (out["saintexpress"]["args"], out["saintexpress"]["cmd_opts"]) == ("-R3 -L50", "-f -x")

    def test_defaults(self):
        assert pfw.generate_saintexpress_config({}) == "-R10 -L100"


class TestCalibratedMzml:
    # FP/FragpipeRun.java:1613,1642: tabMsf.isWriteCalMzml() && tabMsf.getMassCalibration() > 0
    @pytest.mark.parametrize("write,mass,expected", [("true", "2", True), ("true", "0", False), ("false", "2", False)])
    def test_has_calibrated_mzml(self, write, mass, expected):
        p = params(msfragger={"write_calibrated_mzml": write, "calibrate_mass": mass, "table.var-mods": "15.9949,M,true,3"})
        assert pfw.generate_json_output(p, ["DDA"])["msfragger"]["has_calibrated_mzml"] is expected


class TestMinProbJavaRounding:
    # DecimalFormat("#.##") rounds the binary double HALF_EVEN (checked with JDK 17 jshell)
    @pytest.mark.parametrize("text,java", [("0.005", "0.01"), ("0.995", "0.99"), ("0.125", "0.12"), ("0.135", "0.14")])
    def test_rounds_like_java(self, text, java):
        assert pfw.percolator_min_prob(params(percolator={"min-prob": text})) == java
