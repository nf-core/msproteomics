"""bin/plan_to_tool_configs.py and bin/fragpipe_portable.py: the helpers and the small per-tool extractors.

Run: python3 -m unittest discover -s tests/python
Synthetic plans only; the real-plan goldens are in test_plan_to_tool_configs.py.
"""

import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import test_plan_to_tool_configs as base
from test_fragpipe_portable import fp

ptc = base.ptc
pfw = base.load("parse_fragpipe_workflow", base.PARSER)
step, plan = base.step, base.plan


def config(path, text):
    return f'<config path="{path}" size="{len(text)}" sha256="{"0" * 64}" encoding="text"><![CDATA[{text}]]></config>'


def tool(tool_key, *steps, workflow=None):
    configs, _files = ptc.tool_configs(plan(*steps, workflow=workflow))
    return configs[tool_key]


def extract(extractor, cls, name, argvs, configs_xml=""):
    """An extractor called directly on the steps of one class, with the plan holding only those steps."""
    p = plan(*[step(i, cls, name, a, configs_xml) for i, a in enumerate(argvs, 1)])
    return extractor(ptc.Plan(p).steps, ptc.Plan(p))


class PlanHelpersTest(unittest.TestCase):

    def test_strip_options_drops_options_with_values_and_flags(self):
        self.assertEqual(ptc.strip_options(["--a", "1", "--keep", "--raw", "x"], {"--a"}, {"--raw"}), ["--keep", "x"])

    def test_strip_options_value_option_at_the_end_is_dropped(self):
        self.assertEqual(ptc.strip_options(["--keep", "--a"], {"--a"}), ["--keep"])

    def test_same_args_returns_the_common_options(self):
        self.assertEqual(ptc.same_args("t", [["--a"], ["--a"]]), ["--a"])

    def test_same_args_differing_variants_stop(self):
        with self.assertRaisesRegex(ptc.PlanError, "t: per-file steps pass different options"):
            ptc.same_args("t", [["--a"], ["--b"]])

    def test_by_name_selects_steps_of_that_name(self):
        steps = ET.fromstring(f"<plan>{step(1, 'C', 'A', ['x'])}{step(2, 'C', 'B', ['y'])}</plan>").findall("step")
        self.assertEqual([s.get("order") for s in ptc.by_name(steps, "B")], ["2"])

    def test_key_values_skips_comments_and_trailing_comments(self):
        self.assertEqual(ptc.key_values("# c = 1\na = 1 # note\nb=two\nno equals\n"), {"a": "1", "b": "two"})

    def test_argv_without_argv_element_stops(self):
        with self.assertRaisesRegex(ptc.PlanError, "has no <argv>"):
            ptc.argv(ET.fromstring('<step order="3" name="X"/>'))

    def test_tokenizer_puts_paths_under_the_longest_root_at_path_boundaries(self):
        tok = ptc.Tokenizer({"{WORKDIR}": "/run/fp_work", "{DB_DIR}": "/run", "{LCMS_DIR:1}": "/run/raw/p1"})
        self.assertEqual(tok.path("/run/fp_work/a /run/db.fasta /run/raw/p1/x.mzML /run/fp_workx"),
                         "{WORKDIR}/a {DB_DIR}/db.fasta {LCMS_DIR:1}/x.mzML {DB_DIR}/fp_workx")

    def test_plan_config_missing_stops_and_config_named_finds_by_file_name(self):
        p = ptc.Plan(plan(step(1, "CmdCrystalc", "Crystal-C", ["x"], config("{WORKDIR}/sub/a.params", "k = 1"))))
        self.assertEqual(p.config_named("a.params"), [("{WORKDIR}/sub/a.params", "k = 1")])
        with self.assertRaisesRegex(ptc.PlanError, "the plan has no config nope"):
            p.config("nope")

    def test_saved_workflow_unescapes_java_properties(self):
        p = ptc.Plan(plan(workflow="a.path=C\\:\\\\x\\=y\n# c=1\n"))
        self.assertEqual(ptc.saved_workflow(p)["a.path"], "C:\\x=y")

    def test_saved_workflow_undoes_every_escape_properties_store_writes(self):
        # Properties.saveConvert (Properties.store): \\#, \\!, \\t, \\n, \\r, \\f, a leading "\\ " and \\uXXXX
        text = ("a.opts=--x \\#1 \\!2\\tT\\nN\\rR\\fF\n"
                "b.lead=\\  two leading\n"
                "c.uni=\\u00e9t\\u00e9\n"
                "d.cont=first \\\n    second\n")
        props = ptc.saved_workflow(ptc.Plan(plan(workflow=text)))
        self.assertEqual((props["a.opts"], props["b.lead"], props["c.uni"], props["d.cont"]),
                         ("--x #1 !2\tT\nN\rR\fF", "  two leading", "\u00e9t\u00e9", "first second"))

    def test_saved_workflow_reads_as_the_workflow_parser_does(self):
        text = "x\\:y=a\\=b c\nk v\n! bang comment\nu=\\u0041\\t\nlast=1\nlast=2\n"
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "w.workflow"
            f.write_text(text, encoding="utf-8")
            self.assertEqual(ptc.workflow_properties(text), pfw.load_java_properties(str(f)))

    def test_saved_workflow_needs_exactly_one_config(self):
        root = plan()
        root.remove(next(c for c in root.iter("config")))
        with self.assertRaisesRegex(ptc.PlanError, "has 0 saved fragpipe.workflow configs"):
            ptc.saved_workflow(ptc.Plan(root))

    def test_workflow_value_missing_key_stops(self):
        with self.assertRaisesRegex(ptc.PlanError, "has no tmtintegrator.channel_num"):
            ptc.workflow_value(ptc.Plan(plan()), "tmtintegrator.channel_num")

    def test_workflow_value_reads_the_saved_key(self):
        self.assertEqual(ptc.workflow_value(ptc.Plan(plan(workflow="x.y=7\n")), "x.y"), "7")


class SettingsOnlyToolsTest(unittest.TestCase):

    def test_abacus_off_keeps_the_filters_picked_and_adds_reprint_for_saintexpress(self):
        props = {"phi-report.filter": "--sequential --picked --prot 0.01", "database.decoy-tag": "rev_",
                 "phi-report.dont-use-prot-proph-file": "false", "phi-report.prot-level-summary": "true",
                 "phi-report.pep-level-summary": "TRUE"}
        entry = ptc.abacus_off_config(props)
        self.assertEqual(entry["args"], "--picked --razor --reprint --tag rev_ --protein --peptide")
        self.assertFalse(entry["run"])

    def test_abacus_off_without_protein_inference_has_no_protein_summary(self):
        props = {"phi-report.filter": "--sequential", "database.decoy-tag": "rev_", "phi-report.dont-use-prot-proph-file": "true",
                 "phi-report.prot-level-summary": "true", "phi-report.pep-level-summary": "false"}
        self.assertEqual(ptc.abacus_off_config(props)["args"], "--razor --reprint --tag rev_")

    def test_saintexpress_defaults_and_cmd_opts(self):
        self.assertEqual(ptc.saintexpress_config({})["args"], "-R10 -L100")
        entry = ptc.saintexpress_config({"saintexpress.max-replicates": "3", "saintexpress.virtual-controls": "50", "saintexpress.cmd-opts": " -x "})
        self.assertEqual((entry["args"], entry["cmd_opts"]), ("-R3 -L50", "-x"))

    def test_saintexpress_runs_only_when_the_saved_workflow_asks(self):
        self.assertFalse(tool("saintexpress", workflow="saintexpress.run-saint-express=false\n")["run"])
        self.assertTrue(tool("saintexpress", workflow="saintexpress.run-saint-express=true\n")["run"])

    def test_database_decoy_tag_must_match_the_annotation_step(self):
        annotate = step(1, "CmdPhilosopherDbAnnotate", "Annotate", ["philosopher", "database", "--prefix", "decoy_"])
        with self.assertRaisesRegex(ptc.PlanError, "differs from database.decoy-tag=rev_"):
            ptc.tool_configs(plan(annotate))

    def test_database_entry_carries_the_decoy_tag(self):
        self.assertEqual(tool("database"), {"run": True, "args": "--prefix rev_", "config_type": "cli", "decoy_tag": "rev_"})

    def test_a_step_without_extractor_stops(self):
        with self.assertRaisesRegex(ptc.PlanError, "CmdUnknown, X\\) has no extractor yet"):
            ptc.tool_configs(plan(step(1, "CmdUnknown", "X", ["a"])))


class ProphetAndPhilosopherExtractorsTest(unittest.TestCase):

    def test_proteinprophet_drops_output_option_and_the_filelist(self):
        entry = tool("proteinprophet", step(1, "CmdProteinProphet", "ProteinProphet",
                                             ["philosopher", "proteinprophet", "--maxprobability", "--output", "combined", "filelist.txt"]))
        self.assertEqual(entry["args"], "--maxprobability")

    def test_filter_and_report_args(self):
        entry = tool(
            "filter",
            step(1, "CmdPhilosopherFilter", "PhilosopherFilter",
                 ["philosopher", "filter", "--sequential", "--prot", "0.01", "--tag", "rev_", "--pepxml", "d", "--razor", "--protxml", "p.xml"]),
            step(2, "CmdPhilosopherReport", "PhilosopherReport", ["philosopher", "report", "--msstats"]),
        )
        self.assertEqual((entry["args"], entry["report_args"], entry["dont_use_protxml"]), ("--sequential --prot 0.01", "--msstats", False))

    def test_filter_without_protxml_marks_dont_use_protxml(self):
        entry = tool(
            "filter",
            step(1, "CmdPhilosopherFilter", "PhilosopherFilter", ["philosopher", "filter", "--pepxml", "d"]),
            step(2, "CmdPhilosopherReport", "PhilosopherReport", ["philosopher", "report"]),
        )
        self.assertTrue(entry["dont_use_protxml"])

    def test_ptmprophet_args_drop_threads_and_files(self):
        entry = tool("ptmprophet", step(1, "CmdPtmProphet", "PTMProphet", ["PTMProphetParser", "STY:79.966", "MAXTHREADS=1", "a.pepXML", "a.mod.pep.xml"]))
        self.assertEqual(entry["args"], "STY:79.966")

    def test_ptmprophet_argv_not_ending_with_the_files_stops(self):
        with self.assertRaisesRegex(ptc.PlanError, "does not end with MAXTHREADS=1"):
            tool("ptmprophet", step(1, "CmdPtmProphet", "PTMProphet", ["PTMProphetParser", "STY:79.966", "a.pepXML"]))

    def test_iprophet_args_drop_threads_and_pepxmls(self):
        entry = tool("iprophet", step(1, "CmdIprophet", "iProphet",
                                       ["philosopher", "iprophet", "--decoy", "rev_", "--nonsp", "--threads", "4", "a.pep.xml", "b.pepXML"]))
        self.assertEqual(entry["args"], "--decoy rev_ --nonsp")

    def test_iprophet_needs_one_step(self):
        two = [step(i, "CmdIprophet", "iProphet", ["philosopher", "iprophet"]) for i in (1, 2)]
        with self.assertRaisesRegex(ptc.PlanError, "expected one iProphet step, plan has 2"):
            tool("iprophet", *two)

    def test_labelquant_args_and_brand(self):
        entry = tool("labelquant", step(1, "CmdLabelquant", "LabelQuant",
                                         ["philosopher", "labelquant", "--tol", "20", "--brand", "tmt", "--annot", "x.txt", "--raw", "--dir", "d"]))
        self.assertEqual((entry["args"], entry["brand"]), ("--tol 20", "tmt"))

    def test_labelquant_without_brand_stops(self):
        with self.assertRaisesRegex(ptc.PlanError, "has no --brand"):
            tool("labelquant", step(1, "CmdLabelquant", "LabelQuant", ["philosopher", "labelquant", "--tol", "20"]))

    def test_freequant_before_labelquant_does_not_run_label_free(self):
        freequant = step(1, "CmdFreequant", "FreeQuant", ["philosopher", "freequant", "--ptw", "0.4", "--dir", "d"])
        labelquant = step(2, "CmdLabelquant", "LabelQuant", ["philosopher", "labelquant", "--brand", "tmt"])
        configs, _files = ptc.tool_configs(plan(freequant, labelquant))
        self.assertEqual((configs["freequant"]["run"], configs["freequant"]["args"]), (False, "--ptw 0.4"))

    def test_freequant_alone_runs(self):
        entry = tool("freequant", step(1, "CmdFreequant", "FreeQuant", ["philosopher", "freequant", "--ptw", "0.4", "--raw", "--dir", "d"]))
        self.assertEqual((entry["run"], entry["args"]), (True, "--ptw 0.4"))


class SmallToolExtractorsTest(unittest.TestCase):

    def test_appendfile_runs_without_args(self):
        self.assertEqual(tool("appendfile", step(1, "CmdAppendFile", "AppendFile", ["x"]))["args"], "")

    def test_writesubmzml_threshold(self):
        entry = tool("writesubmzml", step(1, "CmdWriteSubMzml", "Write sub mzML",
                                           ["java", "-cp", "x", "a.b.WriteSubMzml", "a.mzML", "psm.tsv", "a_sub.mzML", "0.99", "1"]))
        self.assertEqual(entry["prob_threshold"], "0.99")

    def test_writesubmzml_different_thresholds_stop(self):
        steps = [step(i, "CmdWriteSubMzml", "Write sub mzML", ["java", "a.b.WriteSubMzml", "a.mzML", "p", "s", t, "1"]) for i, t in ((1, "0.9"), (2, "0.8"))]
        with self.assertRaisesRegex(ptc.PlanError, "different probability thresholds"):
            tool("writesubmzml", *steps)

    def test_mbg_runtime_options_become_tokens_and_tools_dir_a_token(self):
        argv = ["java", "-cp", "x", "com.mbg.MBG", "--psm", "g/psm.tsv", "--manifest", "m.fp-manifest", "--glycans", "{FRAGPIPE}/lib/../tools/g.csv", "--numthreads", "4"]
        entry = tool("mbg", step(1, "CmdMBGMatch", "MBG", argv), workflow="workflow.input.data-type.im-ms=TRUE\n")
        self.assertEqual(entry["args"], ["--psm", "@PSM@", "--manifest", "@MANIFEST@", "--glycans", "@TOOLS_DIR@/g.csv", "--numthreads", "@THREADS@"])
        self.assertTrue(entry["im_ms"])

    def test_mbg_step_without_the_class_stops(self):
        with self.assertRaisesRegex(ptc.PlanError, "does not run com.mbg.MBG"):
            tool("mbg", step(1, "CmdMBGMatch", "MBG", ["java", "other"]), workflow="workflow.input.data-type.im-ms=false\n")

    def test_pairscans_args_follow_the_file_and_thread_count(self):
        entry = tool("pairscans", step(1, "CmdPairScans", "PairScans", ["java", "-cp", "x", "org.nesvilab.fragpipe.util.PairScans", "a.mzML", "4", "--min", "1"]))
        self.assertEqual(entry["args"], ["--min", "1"])

    def test_fpop_coadaptr_arg_is_the_last_argv_element(self):
        self.assertEqual(tool("fpop_coadaptr", step(1, "CmdFPOPcoadaptr", "coADAPTr", ["java", "FPOP_PSMwriter", "/w", "offsets.txt"]))["args"], "offsets.txt")

    def test_fpop_tmt_tables_keep_their_report_dirs_in_the_work_dir(self):
        # CmdFpopQuant.java:91-98, :114-145: <wd>/tmt-report/abundance_multi-site_<norm>.tsv first, the
        # unmodified <wd>/tmt-report-unmod/abundance_peptide_<norm>.tsv last, both read in the FPOP task's work dir
        entry = tool("fpop", step(1, "CmdFpopQuant", "FPOP", ["python", "/x/FragPipe_FPOP_Analysis.py", "{WORKDIR}/tmt-report/abundance_multi-site_MD.tsv",
                                                              "1", "ctrl", "fpop", "true", "true", "{WORKDIR}/tmt-report-unmod/abundance_peptide_MD.tsv"]))
        self.assertEqual(entry["args"], ["tmt-report/abundance_multi-site_MD.tsv", "1", "ctrl", "fpop", "true", "true",
                                         "tmt-report-unmod/abundance_peptide_MD.tsv"])

    def test_fpop_label_free_args_follow_the_script_with_the_table_in_the_work_dir(self):
        entry = tool("fpop", step(1, "CmdFpopQuant", "FPOP", ["python", "/x/FragPipe_FPOP_Analysis.py", "{WORKDIR}/combined_modified_peptide.tsv",
                                                              "1", "ctrl", "fpop", "true", "false"]))
        self.assertEqual(entry["args"], ["combined_modified_peptide.tsv", "1", "ctrl", "fpop", "true", "false"])

    def test_metaproteomics_options_and_trailing_opts(self):
        opts = {"qvalue": "0.01", "deltaHyperscore": "3", "minPeptCntPerProt": "2", "minUniqPeptCntPerProt": "1", "minUniqPeptCnt": "1", "hostName": "human", "iterations": "5"}
        argv = ["java", "-cp", "x", "DbOptimizer"] + [x for k, v in opts.items() for x in (f"--{k}", v)] + ["--extra", "1"]
        args = tool("metaproteomics", step(1, "CmdMetaproteomics", "FP-Meta", argv))["args"]
        self.assertEqual((args["qvalue"], args["hostName"], args["cmdLineOpts"]), ("0.01", "human", "--extra 1"))

    def test_crystalc_args_drop_comments_and_runtime_keys(self):
        text = "# header\nthread = 4\nfasta = /x.fasta\nprecursorIsolationWindow = 0.7 # comment\n\nmaxCharge = 3\n"
        entry = extract(ptc.crystalc_config, "CmdCrystalc", "Crystal-C", [["java"]], config("{WORKDIR}/crystalc-0-a.params", text))
        self.assertEqual(entry["args"], "precursorIsolationWindow = 0.7\nmaxCharge = 3")

    def test_crystalc_needs_one_params_file_per_step(self):
        with self.assertRaisesRegex(ptc.PlanError, "has 0 crystalc params files"):
            extract(ptc.crystalc_config, "CmdCrystalc", "Crystal-C", [["java"]])

    def test_diaumpire_threads_token_and_q_tiers(self):
        umpire = step(1, "CmdUmpireSe", "DIA-Umpire SE", ["java"], config("{WORKDIR}/umpire-se.params", "Thread = 8\nRPmax = 25"))
        search = step(2, "CmdMsfragger", "MSFragger", ["java", "a_Q1.mzML", "a_Q3.mzML", "a.mzML"])
        p = ptc.Plan(plan(umpire, search))
        entry = ptc.diaumpire_config([s for s in p.steps if s.get("class") == "CmdUmpireSe"], p)
        self.assertEqual((entry["args"], entry["q_tiers"]), ("Thread = THREADS\nRPmax = 25", ["Q1", "Q3"]))

    def test_diaumpire_keeps_the_files_final_newlines(self):
        # FragPipe println()s every template line, so umpire-se.params ends with an empty line (PropertiesUtils.java:291-321)
        umpire = step(1, "CmdUmpireSe", "DIA-Umpire SE", ["java"], config("{WORKDIR}/umpire-se.params", "Thread = 8\nRPmax = 25\n\n"))
        p = ptc.Plan(plan(umpire))
        entry = ptc.diaumpire_config([s for s in p.steps if s.get("class") == "CmdUmpireSe"], p)
        self.assertEqual(entry["args"], "Thread = THREADS\nRPmax = 25\n\n")

    def test_diaumpire_needs_one_params_file(self):
        with self.assertRaisesRegex(ptc.PlanError, "expected one umpire-se.params, plan has 0"):
            extract(ptc.diaumpire_config, "CmdUmpireSe", "DIA-Umpire SE", [["java"]])

    def test_opair_runtime_options_and_calibrated_spectra(self):
        argv = ["dotnet", "CMD.dll", "-t 4", "-r a", "-s b", "-o c", "-g {FRAGPIPE}/lib/../tools/glycan.csv"]
        entry = extract(ptc.opair_config, "CmdOPair", "O-Pair", [argv], config("{WORKDIR}/filelist_opair.txt", "a_calibrated.mzML\n"))
        self.assertEqual(entry["args"], ["-t @THREADS@", "-r @FILELIST@", "-s @PSM@", "-o @OUTDIR@", "-g @TOOLS_DIR@/glycan.csv"])
        self.assertEqual(entry["spectra"], "calibrated")

    def test_opair_uncalibrated_spectra_and_missing_filelist(self):
        entry = extract(ptc.opair_config, "CmdOPair", "O-Pair", [["dotnet", "CMD.dll"]], config("{WORKDIR}/filelist_opair.txt", "a.mzML\n"))
        self.assertEqual(entry["spectra"], "uncalibrated_or_input")
        with self.assertRaisesRegex(ptc.PlanError, "no filelist_opair.txt"):
            extract(ptc.opair_config, "CmdOPair", "O-Pair", [["dotnet", "CMD.dll"]])


class PortableHelpersTest(unittest.TestCase):

    def test_boundary_matches_whole_paths_and_directory_prefixes_only(self):
        pattern = fp.boundary("/w/a")
        self.assertTrue(pattern.search("x=/w/a"))
        self.assertTrue(pattern.search("/w/a/b.txt"))
        self.assertFalse(pattern.search("/w/ab"))


if __name__ == "__main__":
    unittest.main()
