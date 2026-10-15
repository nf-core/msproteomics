"""Quantification options of parse_fragpipe_workflow.py against FragPipe 24.0.

Run: make test-python
Expected values come from the FragPipe 24.0 sources cited in each test and from the golden
FragPipe run (TMT-labelcheck-304_perfile, IonQuant step argv).
"""

import importlib.util
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "modules/local/parsefragpipeworkflow/resources/usr/bin/parse_fragpipe_workflow.py"
# FragPipe 24.0 templates the module stages (the script takes them as arguments)
TEMPLATES = Path(__file__).resolve().parents[2] / "modules" / "local" / "parsefragpipeworkflow" / "assets"
TEMPLATE_ARGS = ["--fragger-template", str(TEMPLATES / "fragger_closed.params"), "--diaumpire-template", str(TEMPLATES / "diaumpire_se.params")]
FIXTURES = REPO / "tests/fixtures/fragpipe_workflows"

spec = importlib.util.spec_from_file_location("parse_fragpipe_workflow", SCRIPT)
pfw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pfw)
pfw.FRAGGER_TEMPLATE = TEMPLATES / "fragger_closed.params"
pfw.DIAUMPIRE_TEMPLATE = TEMPLATES / "diaumpire_se.params"


def parse_with_threads(path: str) -> dict:
    """parse_workflow_file plus the run's thread count: the pipeline passes --threads (workflow.threads), which
    the parser requires whenever MSBooster runs (its RT loess output depends on it)."""
    params = pfw.parse_workflow_file(path)
    params.setdefault("workflow", {})["threads"] = "4"
    return params


def workflow(name, **overrides):
    """Parsed bundled workflow with key overrides ("prefix.key" -> value)."""
    params = parse_with_threads(str(FIXTURES / f"{name}.workflow"))
    for key, value in overrides.items():
        prefix, sub = key.split(".", 1)
        params.setdefault(prefix, {})[sub] = value
    return params


# Golden FragPipe 24.0 IonQuant argv (TMT-labelcheck-304 workflow) between --threads and --multidir
GOLDEN_LFQ_ARGS = (
    "--perform-ms1quant 1 --perform-isoquant 0 --isotol 20.0 --isolevel 2 --isotype tmt10 "
    "--ionmobility 0 --site-reports 1 --msstats 1 --minexps 1 --mbr 0 --maxlfq 0 --requantify 1 "
    "--mztol 10 --imtol 0.05 --rttol 0.4 --mbrmincorr 0 --mbrrttol 1 --mbrimtol 0.05 --mbrtoprun 10 "
    "--ionfdr 0.01 --proteinfdr 1 --peptidefdr 1 --normalization 0 --minisotopes 2 --intensitymode 2 "
    "--minscans 3 --writeindex 0 --tp 0 --minfreq 0 --minions 1 --locprob 0.75 --uniqueness 0"
)


class TestIonQuantLfq(unittest.TestCase):
    def test_lfq_args_match_golden_fragpipe_argv(self):
        params = parse_with_threads(str(REPO / "assets/TMT-labelcheck-304.workflow"))
        self.assertEqual(pfw.generate_ionquant_config(params["ionquant"], params), GOLDEN_LFQ_ARGS)

    def test_msstats_is_always_one_for_lfq(self):
        # FragpipeRun.java:1888 passes generateMSstatsInput=true; diann.generate-msstats is DIA-NN only
        params = workflow("LFQ-MBR", **{"diann.generate-msstats": "false"})
        self.assertIn("--msstats 1", pfw.generate_ionquant_config(params["ionquant"], params))

    def test_minexps_and_ionmobility_ignore_workflow_keys(self):
        # CmdIonquant.java:197-198, 216-219
        params = workflow("LFQ-MBR", **{"ionquant.minexps": "3", "ionquant.ionmobility": "1"})
        args = pfw.generate_ionquant_config(params["ionquant"], params).split()
        self.assertEqual(args[args.index("--minexps") + 1], "1")
        self.assertEqual(args.count("--minexps"), 1)
        self.assertEqual(args[args.index("--ionmobility") + 1], "0")
        self.assertEqual(args.count("--ionmobility"), 1)

    def test_ionmobility_from_im_ms_button(self):
        params = workflow("LFQ-MBR", **{"workflow.input.data-type.im-ms": "true"})
        self.assertIn("--ionmobility 1", pfw.generate_ionquant_config(params["ionquant"], params))

    def test_labels_only_with_use_labeling(self):
        params = workflow("LFQ-MBR", **{"ionquant.light": "K0", "ionquant.use-labeling": "false"})
        self.assertNotIn("--light", pfw.generate_ionquant_config(params["ionquant"], params))
        params["ionquant"]["use-labeling"] = "true"
        self.assertIn("--light K0", pfw.generate_ionquant_config(params["ionquant"], params))


class TestIonQuantPanelDefaults(unittest.TestCase):
    def test_missing_keys_get_widget_defaults(self):
        # QuantPanelLabelfree.java:274-354 initial values
        params = workflow("LFQ-MBR")
        for key in ("mbr", "mztol", "rttol", "minscans", "locprob", "uniqueness"):
            params["ionquant"].pop(key)
        args = pfw.generate_ionquant_config(params["ionquant"], params)
        for opt in ("--mbr 1", "--mztol 10", "--rttol 0.4", "--minscans 3", "--locprob 0.75", "--uniqueness 0"):
            self.assertIn(opt, args)

    def test_values_are_formatted_like_the_spinners(self):
        params = workflow("LFQ-MBR", **{"ionquant.rttol": "1", "ionquant.mztol": "10.0", "ionquant.imtol": "0.1", "ionquant.mbr": "true"})
        args = pfw.generate_ionquant_config(params["ionquant"], params)
        for opt in ("--rttol 1.0", "--mztol 10 ", "--imtol 0.10", "--mbr 1"):
            self.assertIn(opt, args)

    def test_invalid_combo_value_fails(self):
        params = workflow("LFQ-MBR", **{"ionquant.intensitymode": "peak"})
        with self.assertRaises(ValueError):
            pfw.generate_ionquant_config(params["ionquant"], params)


class TestModMasses(unittest.TestCase):
    def test_golden_modmasses_file(self):
        params = parse_with_threads(str(REPO / "assets/TMT-labelcheck-304.workflow"))
        # golden work/modmasses_ionquant.txt
        self.assertEqual(pfw._parse_mod_masses(params), ["15.9949", "42.0106", "57.02146", "304.2071"])

    def test_java_float_to_string_of_tmt6_mass(self):
        params = workflow("LFQ-MBR", **{"msfragger.table.var-mods": "229.162932,K,true,1", "msfragger.table.fix-mods": ""})
        self.assertEqual(pfw._parse_mod_masses(params), ["229.16293"])

    def test_rows_without_sites_are_skipped(self):
        # TabMsfragger.java:1709: mod.sites != null && !mod.sites.trim().isEmpty()
        params = workflow("LFQ-MBR", **{"msfragger.table.var-mods": "79.96633, ,true,3; 15.9949,M,false,3", "msfragger.table.fix-mods": ""})
        self.assertEqual(pfw._parse_mod_masses(params), [])

    def test_no_modlist_when_opair_runs(self):
        params = workflow("LFQ-MBR", **{"opair.run-opair": "true"})
        self.assertIsNone(pfw._parse_mod_masses(params))

    def test_tmt_passes_get_modmasses_without_lfq_ionquant(self):
        params = workflow("TMT10", **{"ionquant.run-ionquant": "false"})
        out = pfw.generate_json_output(params, ["DDA"])
        self.assertEqual(out["ionquant_ms1"]["modmasses"], pfw._parse_mod_masses(params))
        self.assertTrue(out["ionquant_ms1"]["modmasses"])


class TestTmtIonQuantPasses(unittest.TestCase):
    def test_ms1_pass_fixed_options(self):
        # CmdIonquant.java:178-254 with uiCompsRepresentation == null
        out = pfw.generate_json_output(workflow("TMT10"), ["DDA"])
        self.assertEqual(
            out["ionquant_ms1"]["args"],
            "--perform-ms1quant 1 --perform-isoquant 0 --isotol 20.0 --isolevel 2 --isotype tmt10 "
            "--ionmobility 0 --site-reports 0 --msstats 0 --minexps 1 --mbr 0 --maxlfq 0 --requantify 0 "
            "--mztol 10 --imtol 0.05 --rttol 1 --normalization 0 --minisotopes 1 --minscans 1 "
            "--writeindex 0 --tp 0 --minfreq 0 --minions 1 --locprob 0 --uniqueness 0",
        )

    def test_isobaric_pass_label_and_level(self):
        out = pfw.generate_json_output(workflow("TMT10", **{"tmtintegrator.quant_level": "3", "tmtintegrator.tolerance": "15"}), ["DDA"])
        args = out["ionquant_isobaric"]["args"]
        self.assertTrue(args.startswith("--perform-ms1quant 0 --perform-isoquant 1 --isotol 15.0 --isolevel 3 --isotype TMT-10 "))

    def test_no_ionquant_passes_for_philosopher(self):
        out = pfw.generate_json_output(workflow("TMT10", **{"tmtintegrator.extraction_tool": "Philosopher", "freequant.run-freequant": "false"}), ["DDA"])
        self.assertNotIn("ionquant_ms1", out)
        self.assertNotIn("ionquant_isobaric", out)

    def test_skip_extraction_fails(self):
        # FragpipeRun.java:1913-1916 (Report runs in TMT10)
        with self.assertRaises(ValueError):
            pfw.generate_json_output(workflow("TMT10", **{"tmtintegrator.extraction_tool": "Skip extraction. Run TMT-Integrator only"}), ["DDA"])

    def test_philosopher_msstats_needs_philosopher_extraction(self):
        # FragpipeRun.java:1768-1773
        with self.assertRaises(ValueError):
            pfw.generate_json_output(workflow("TMT10", **{"tmtintegrator.philosopher-msstats": "true"}), ["DDA"])
        out = pfw.generate_json_output(workflow("TMT10", **{
            "tmtintegrator.philosopher-msstats": "true", "tmtintegrator.extraction_tool": "Philosopher", "freequant.run-freequant": "false"}), ["DDA"])
        self.assertTrue(out["tmtintegrator"]["philosopher_msstats"])

    def test_min_site_prob_needs_mod_tag(self):
        # CmdTmtIntegrator.java:245-253
        with self.assertRaises(ValueError):
            pfw.generate_json_output(workflow("TMT10", **{"tmtintegrator.min_site_prob": "0.75", "tmtintegrator.mod_tag": ""}), ["DDA"])
        pfw.generate_json_output(workflow("TMT10", **{"tmtintegrator.min_site_prob": "0.75", "tmtintegrator.mod_tag": "S(79.9663)"}), ["DDA"])

    def test_unknown_extraction_tool_fails(self):
        with self.assertRaises(ValueError):
            pfw.generate_json_output(workflow("TMT10", **{"tmtintegrator.extraction_tool": "1"}), ["DDA"])


class TestPhilosopherExtraction(unittest.TestCase):
    def philosopher(self, **overrides):
        return workflow("TMT10", **{"tmtintegrator.extraction_tool": "Philosopher", "freequant.run-freequant": "false", **overrides})

    def test_labelquant_runs_for_philosopher_text(self):
        out = pfw.generate_json_output(self.philosopher(), ["DDA"])
        self.assertTrue(out["labelquant"]["run"])
        self.assertEqual(out["tmtintegrator"]["extraction_tool"], 1)

    def test_labelquant_options(self):
        # CmdLabelquant.java:98-116 (String.format("%.3f")), plex = reagent count
        out = pfw.generate_json_output(self.philosopher(**{"tmtintegrator.min_pep_prob": "0.9", "tmtintegrator.min_percent": "0.05"}), ["DDA"])
        self.assertEqual(out["labelquant"]["args"], "--tol 20 --level 2 --minprob 0.900 --purity 0.500 --removelow 0.050 --plex 10")
        self.assertEqual(out["labelquant"]["brand"], "tmt")

    def test_iodotmt_brand_is_label_type(self):
        # QuantLabel.java:59 (TYPE_TMT), CmdLabelquant.java:144-145
        out = pfw.generate_json_output(self.philosopher(**{"tmtintegrator.channel_num": "iodoTMT-6"}), ["DDA"])
        self.assertEqual(out["labelquant"]["brand"], "tmt")
        self.assertIn("--plex 6", out["labelquant"]["args"])

    def test_zoom_hr_rejected(self):
        with self.assertRaises(ValueError):
            pfw.generate_json_output(self.philosopher(**{"tmtintegrator.quant_level": "4"}), ["DDA"])

    def test_freequant_options_for_tmt_freequant(self):
        # FragpipeRun.java:2049-2055 with QuantPanelLabelfree.getFreequantOptsAsText()
        out = pfw.generate_json_output(self.philosopher(**{"freequant.rt-tol": "0.4", "freequant.mz-tol": "10"}), ["DDA"])
        self.assertFalse(out["freequant"]["run"])
        self.assertEqual(out["freequant"]["args"], "--ptw 0.4 --tol 10")

    def test_freequant_and_philosopher_tmt_conflict(self):
        # FreeQuant selected: its twin IonQuant (true in the TMT10 workflow) must not also be true
        params = self.philosopher(**{"freequant.run-freequant": "true", "ionquant.run-ionquant": "false",
                                     "quantitation.run-label-free-quant": "true"})
        if pfw.get_run_flag(params, "freequant"):
            with self.assertRaises(ValueError):
                pfw.generate_json_output(params, ["DDA"])


class TestTmtIntegratorConfig(unittest.TestCase):
    def config(self, **overrides):
        params = workflow("TMT10", **overrides)
        return dict(
            line.strip().split(": ", 1)
            for line in pfw.generate_tmtintegrator_config(params["tmtintegrator"]).splitlines()[1:]
        )

    def test_every_line_ends_with_newline_like_tmticonfig_write(self):
        # TmtiConfig.write (TmtiConfig.java:332-338) ends each line, the last one included, with "\n"
        text = pfw.generate_tmtintegrator_config(workflow("TMT10")["tmtintegrator"])
        self.assertTrue(text.endswith("use_glycan_composition: false\n"))
        self.assertFalse(text.endswith("\n\n"))

    def test_sorted_keys_without_path_or_memory(self):
        params = workflow("TMT10")
        lines = pfw.generate_tmtintegrator_config(params["tmtintegrator"]).splitlines()
        self.assertEqual(lines[0], "tmtintegrator:")
        keys = [line.strip().split(":")[0] for line in lines[1:]]
        self.assertEqual(keys, sorted(keys))
        self.assertNotIn("path", keys)
        self.assertNotIn("memory", keys)
        self.assertTrue(all(line.startswith("  ") and not line.strip().startswith("#") for line in lines[1:]))

    def test_channel_num_and_label_masses(self):
        conf = self.config()
        self.assertEqual((conf["channel_num"], conf["label_masses"], conf["is_tmt_35"]), ("10", "229.16293", "false"))

    def test_ref_tag_is_sanitised(self):
        # TmtiPanel.java:1203-1204, QuantLabelAnnotation.disallowedPattern
        self.assertEqual(self.config(**{"tmtintegrator.ref_tag": "Br idge!"})["ref_tag"], "Bridge")

    def test_output_placeholder(self):
        self.assertEqual(self.config()["output"], pfw.TMTI_OUTPUT_PLACEHOLDER)

    def test_missing_label_is_the_combo_default(self):
        params = parse_with_threads(str(REPO / "assets/TMT-labelcheck-304.workflow"))
        self.assertIn("  channel_num: 1", pfw.generate_tmtintegrator_config(params["tmtintegrator"]).splitlines())

    def test_unknown_label_fails(self):
        with self.assertRaises(ValueError):
            self.config(**{"tmtintegrator.channel_num": "TMT-99"})


if __name__ == "__main__":
    unittest.main()


class TestTmtIntegratorConfigText(unittest.TestCase):
    """tool_configs JSON carries the TMT-Integrator config texts (TMTINTEGRATOR reads no file on the head node):
    tmtintegrator.args and, for FPOP TMT, tmtintegrator_unmod.args equal the files written beside them."""

    def test_first_run_args_equal_tmtintegrator_yml(self):
        params = workflow("TMT10")
        result = pfw.generate_json_output(params, ["DDA"])
        self.assertTrue(result["tmtintegrator"]["run"])
        self.assertEqual(result["tmtintegrator"]["args"], pfw.generate_tmtintegrator_config(params["tmtintegrator"]))

    def test_fpop_tmt_unmod_args_equal_tmtintegrator_unmod_yml(self):
        params = workflow("TMT10-fpop")
        result = pfw.generate_json_output(params, ["DDA"])
        self.assertTrue(result["tmtintegrator_unmod"]["run"])
        self.assertEqual(result["tmtintegrator_unmod"]["args"], pfw.generate_tmtintegrator_config(params["tmtintegrator"], True))
        self.assertIn("  mod_tag: none\n", result["tmtintegrator_unmod"]["args"])

    def test_no_unmod_text_without_fpop_tmt(self):
        self.assertEqual(pfw.generate_json_output(workflow("TMT10"), ["DDA"])["tmtintegrator_unmod"]["args"], "")
