"""bin/plan_to_tool_configs.py: FragPipe's --write-plan work dir to tool_configs.json.

Run: python3 -m unittest discover -s tests/python
Unit tests use small synthetic plans (the plan load_workdir builds) and small synthetic work dirs. The parity
test compares the planner with parse_fragpipe_workflow.py on the goldens: each work dir was written by
FRAGPIPEPLAN (modules/local/fragpipeplan/tests/golden*.nf.test) from the golden run's saved workflow and manifest.
"""

import importlib.util
import json
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PLANS = REPO / "tests/fixtures/fragpipe_plans"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ptc = load("plan_to_tool_configs", REPO / "bin/plan_to_tool_configs.py")
PARSER = REPO / "modules/local/parsefragpipeworkflow/resources/usr/bin/parse_fragpipe_workflow.py"
# FragPipe 24.0 templates the PARSEFRAGPIPEWORKFLOW module stages (the script takes them as arguments)
PARSER_TEMPLATES = ["--fragger-template", str(REPO / "modules/local/parsefragpipeworkflow/assets/fragger_closed.params"),
                    "--diaumpire-template", str(REPO / "modules/local/parsefragpipeworkflow/assets/diaumpire_se.params")]

ROOTS = '<roots><root name="{WORKDIR}" path="/w"/><root name="{FRAGPIPE}" path="/f"/></roots>'


def step(order, cls, name, args, configs=""):
    argv = "".join(f"<arg>{a}</arg>" for a in args)
    return (f'<step order="{order}" tool="t" class="{cls}" name="{name}" gate="SEQUENTIAL EXECUTION">'
            f"<argv>{argv}</argv><inputs/>{configs}</step>")


def config(path, text):
    return f'<config path="{path}" size="1" sha256="{"0" * 64}" encoding="text"><![CDATA[{text}]]></config>'


def saved_workflow(text="workflow.threads=4\n"):
    """The run's fragpipe.workflow FragPipe saves into its work dir: a top-level plan config."""
    return ('<config path="{WORKDIR}/fragpipe.workflow" size="1" sha256="' + "0" * 64 + '" encoding="text">'
            f"<![CDATA[{text}]]></config>")


# keys every saved fragpipe.workflow has (FragPipe saves each UI key) that the planner reads for tools without steps
SAVED_KEYS = ("database.decoy-tag=rev_\nphi-report.filter=--sequential --picked --prot 0.01\nphi-report.dont-use-prot-proph-file=false\n"
              "phi-report.prot-level-summary=true\nphi-report.pep-level-summary=false\n")


def plan(*steps, extra="", workflow=None):
    return ET.fromstring(f'<plan>{ROOTS}{extra}{"".join(steps)}{saved_workflow(SAVED_KEYS + (workflow or "workflow.threads=4\n"))}</plan>')


PERCOLATOR = ["percolator", "--only-psms", "--num-threads", "8", "--results-psms", "t.tsv",
              "--decoy-results-psms", "d.tsv", "--protein-decoy-pattern", "rev_"]
CONVERT = ["java", "-cp", "lib/*", "PercolatorOutputToPepXML", "a.pin", "a", "t.tsv", "d.tsv", "interact-a", "DDA", "0.5", "a.mzML"]


class PlanToToolConfigsTest(unittest.TestCase):

    def test_tools_without_steps_do_not_run(self):
        configs, _files = ptc.tool_configs(plan())
        self.assertEqual(set(configs), set(ptc.CONFIG_TYPES))
        # the database entry always runs: it carries the run's decoy tag (database.decoy-tag), as in the parser
        self.assertEqual([t for t, c in configs.items() if c["run"]], ["database"])

    def test_percolator_args_drop_what_the_module_builds(self):
        configs, _files = ptc.tool_configs(plan(
            step(1, "CmdPercolator", "Percolator", PERCOLATOR + ["a.pin"]),
            step(2, "CmdPercolator", "Percolator: Convert to pepxml", CONVERT),
        ))
        self.assertEqual(configs["percolator"]["args"], "--only-psms")
        self.assertEqual(configs["percolator"]["min_prob"], "0.5")
        self.assertTrue(configs["percolator"]["keep_tsv_files"], "no 'Percolator delete temp' step")

    def test_percolator_write_sub_mzml_weights_are_the_modules(self):
        # first pass: --weights <stem>_percolator.weights per file; second pass: --init-weights <w> --static --override
        configs, _files = ptc.tool_configs(plan(
            step(1, "CmdPercolator", "Percolator", PERCOLATOR + ["--weights", "a_percolator.weights", "a.pin"]),
            step(2, "CmdPercolator", "Percolator", PERCOLATOR + ["--init-weights", "/d/a_percolator.weights", "--static", "--override", "b.pin"]),
            step(3, "CmdPercolator", "Percolator: Convert to pepxml", CONVERT),
        ))
        self.assertEqual(configs["percolator"]["args"], "--only-psms")

    def test_per_file_steps_with_different_options_stop(self):
        with self.assertRaisesRegex(ptc.PlanError, "per-file steps pass different options"):
            ptc.tool_configs(plan(
                step(1, "CmdPercolator", "Percolator", PERCOLATOR + ["a.pin"]),
                step(2, "CmdPercolator", "Percolator", PERCOLATOR + ["--trainFDR", "0.05", "b.pin"]),
                step(3, "CmdPercolator", "Percolator: Convert to pepxml", CONVERT),
            ))

    def test_step_without_extractor_stops(self):
        with self.assertRaisesRegex(ptc.PlanError, "has no extractor yet"):
            ptc.tool_configs(plan(step(1, "CmdSkyline", "Skyline", ["skyline"])))

    def test_export_matched_fragments_runs(self):
        emf = step(1, "CmdExportMatchedFragments", "ExportMatchedFragments", ["java", "-cp", "cp", "FragmentsExportMainClass", "/w", "4", "r"])
        configs, _files = ptc.tool_configs(plan(emf))
        self.assertEqual(configs["exportmatchedfragments"], {"run": True, "args": "", "config_type": "cli"})

    def test_write_sub_mzml_threshold_comes_from_the_step_argv(self):
        def wsm(o, f):
            return step(o, "CmdWriteSubMzml", "WriteSubMzml", ["java", "-Xmx8G", "-cp", "cp", "org.nesvilab.fragpipe.util.WriteSubMzml",
                                                               f"/d/{f}.mzML", "/w/g/psm.tsv", f"/w/{f}_sub.mzML", "0.7", "1"])
        configs, _files = ptc.tool_configs(plan(wsm(1, "a"), wsm(2, "b")))
        self.assertEqual(configs["writesubmzml"], {"run": True, "args": "", "config_type": "cli", "prob_threshold": "0.7"})

    def test_saintexpress_settings_come_from_the_saved_workflow(self):
        configs, _files = ptc.tool_configs(plan(workflow="saintexpress.run-saint-express=true\nsaintexpress.max-replicates=3\n"
                                                         "saintexpress.virtual-controls=100\nsaintexpress.cmd-opts=\n"))
        self.assertEqual(configs["saintexpress"], {"run": True, "args": "-R3 -L100", "config_type": "cli", "cmd_opts": ""})

    def test_saintexpress_off_in_the_saved_workflow_does_not_run(self):
        configs, _files = ptc.tool_configs(plan(workflow="saintexpress.run-saint-express=false\n"))
        self.assertFalse(configs["saintexpress"]["run"])

    def test_split_database_slices_come_from_the_split_script_argv(self):
        params = ('<config path="/w/fragger.params" size="1" sha256="' + "0" * 64 + '" encoding="text">'
                  "<![CDATA[data_type = 0\ncalibrate_mass = 2\nwrite_calibrated_mzml = false\n]]></config>")
        configs, _files = ptc.tool_configs(plan(step(1, "CmdMsfragger", "MSFragger",
                                             ["python3", "/f/msfragger_pep_split.py", "2", "java -jar", "m.jar", "/w/fragger.params"], params)))
        self.assertEqual(configs["msfragger"]["num_slices"], 2)
        self.assertEqual(configs["msfragger"]["calibrate_mass"], {"DDA": "2"})
        self.assertEqual(configs["msfragger"]["args"], "data_type = 0\ncalibrate_mass = 2\nwrite_calibrated_mzml = false\n")

    def test_msbooster_spectral_library_is_staged_through_a_token(self):
        text = "useDetect = false\nfragger = /w/fragger.params\nspectraPredFile = /lib/pred.dlib\nRTPredFile = /lib/pred.dlib\nIMPredFile = /lib/pred.dlib\n"
        configs, _files = ptc.tool_configs(plan(step(1, "CmdMSBooster", "MSBooster", ["java", "-jar", "msbooster.jar"],
                                                     config("/w/msbooster_params.txt", text))))
        self.assertEqual(configs["msbooster"]["spectral_library"], "/lib/pred.dlib")
        self.assertEqual(configs["msbooster"]["args"], "useDetect = false\nspectraPredFile = {SPECTRAL_LIBRARY}\n"
                                                       "RTPredFile = {SPECTRAL_LIBRARY}\nIMPredFile = {SPECTRAL_LIBRARY}")

    def test_ptmshepherd_custom_annotation_file_is_staged_through_a_token(self):
        def shepherd(value):
            text = f"threads = 8\nannotation_file = {value}\n"
            return ptc.tool_configs(plan(step(1, "CmdPtmshepherd", "PTMShepherd", ["java", "-jar", "ptmshepherd.jar", "/w/shepherd.config"],
                                              config("/w/shepherd.config", text))))[0]["ptmshepherd"]
        custom = shepherd("/anno/my_mods.txt")
        self.assertEqual((custom["annotation_file"], custom["args"]), ("/anno/my_mods.txt", "threads = @THREADS@\nannotation_file = @ANNOTATION_FILE@"))
        for keyword in ("unimod", "common", "glyco", ""):
            keyword_entry = shepherd(keyword)
            self.assertEqual((keyword_entry["annotation_file"], keyword_entry["args"]), ("", f"threads = @THREADS@\nannotation_file = {keyword}"))

    def test_opair_user_oxonium_file_is_staged_through_a_token(self):
        def opair(oxo):
            return ptc.tool_configs(plan(step(1, "CmdOPair", "OPair", ["dotnet", "{FRAGPIPE}/tools/opair/CMD.dll", "-b 20.0", f"-f {oxo}", "-m 0.05",
                                                                     "-t 8", "-r /w/g/filelist_opair.txt", "-s /w/g/psm.tsv", "-o /w/g"],
                                              config("/w/g/filelist_opair.txt", "a_calibrated.mzML\n"))))[0]["opair"]
        user = opair("/oxo/my_oxonium.txt")
        self.assertEqual((user["oxonium_file"], user["args"]), ("/oxo/my_oxonium.txt", ["-b 20.0", "-f @OXONIUM_FILE@", "-m 0.05", "-t @THREADS@", "-r @FILELIST@", "-s @PSM@", "-o @OUTDIR@"]))
        default = opair("{FRAGPIPE}/lib/../tools/opair/Glycan_Mods/oxonium_filter.txt")
        self.assertEqual(default["oxonium_file"], "")
        self.assertIn("-f @TOOLS_DIR@/opair/Glycan_Mods/oxonium_filter.txt", default["args"])

    def test_diann_cmd_opts_is_its_own_field(self):
        diann = ["diann", "--lib", "library.tsv", "--threads", "8", "--qvalue", "0.01", "--matrices", "--no-prot-inf",
                 "--smart-profiling", "--peak-center", "--no-ifs-removal", "--report-lib-info", "--cfg", "/w/filelist_diann.txt--"]
        def run(cmd_opts):
            return ptc.tool_configs(plan(step(1, "CmdDiann", "DIA-NN: run DIA-NN", diann[:-2] + ([cmd_opts] if cmd_opts else []) + diann[-2:]),
                                         workflow=f"diann.cmd-opts={cmd_opts}\n"))[0]["diann"]
        entry = run("--mass-acc 15")
        self.assertEqual((entry["cmd_opts"], entry["args_post"]), ("--mass-acc 15", "--smart-profiling --peak-center --no-ifs-removal --report-lib-info"))
        self.assertNotIn("predicted_library", entry)
        self.assertEqual(run("")["cmd_opts"], "")
        with self.assertRaisesRegex(ptc.PlanError, "not the argv element before --cfg"):
            ptc.tool_configs(plan(step(1, "CmdDiann", "DIA-NN: run DIA-NN", diann), workflow="diann.cmd-opts=--mass-acc 15\n"))



JOB_HEADER = "# Job name (optional)\tWorkflow file path\tManifest file path\tOutput directory\tTools folder path (optional)\tFasta path (optional)\tRAM (optional)\tThreads (optional)\n"


def dry_run_log(*steps, lcms="/run/raw/a.mzML"):
    """A FragPipe dry-run log (FragpipeRun.java:368, :412, :477): steps are (header, printed command) pairs."""
    body = "".join(f"{h}\n{c}\n" for h, c in steps)
    return (f"LCMS files:\n  Experiment/Group: a\n  - {lcms}\tDDA\n\n\n{len(steps)} commands to execute:\n{body}"
            "~~~~~~~~~~~~~~~~~~~~~~\n\nIt's a dry-run, not running the commands.\n")


def workdir(tmp, *steps, files=None, tools="/fp/fragpipe-24.0/tools", workflow=""):
    """A --write-plan work dir of the run /run/fp_work: run log, fragpipe.job, fragpipe.workflow, files."""
    wd = Path(tmp) / "fp_work"
    wd.mkdir()
    (wd / "log_2026-10-12_10-00-00.txt").write_text(dry_run_log(*steps))
    (wd / "fragpipe.job").write_text(JOB_HEADER + f"runtime-save\t/run/fp_work/fragpipe-workflow.workflow\t/run/fp_work/fragpipe-files.fp-manifest\t/run/fp_work\t{tools}\t/run/db/db.fasta\t8\t4\n")
    (wd / "fragpipe.workflow").write_text(SAVED_KEYS + "workflow.threads=4\n" + workflow)
    (wd / "fragpipe-files.fp-manifest").write_text("/run/raw/a.mzML\ta\t\tDDA")
    for name, text in (files or {}).items():
        (wd / name).parent.mkdir(parents=True, exist_ok=True)
        (wd / name).write_text(text)
    return wd


MSBOOSTER_STEP = ("MSBooster [Work dir: /run/fp_work]",
                  "java -Xmx8G -cp /fp/fragpipe-24.0/tools/msbooster-1.4.17.jar Features --paramsList /run/fp_work/msbooster_params.txt")


class LoadWorkdirTest(unittest.TestCase):
    """load_workdir: the plan of a --write-plan work dir."""

    def test_steps_roots_and_configs_come_from_the_log_and_the_work_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = ptc.load_workdir(workdir(tmp, MSBOOSTER_STEP, files={"msbooster_params.txt": "useRT = true\n", "g/filelist.txt": "x\n"}))
        roots = {r.get("name"): r.get("path") for r in root.find("roots")}
        self.assertEqual(roots, {"{WORKDIR}": "/run/fp_work", "{LCMS_DIR}": "/run/raw", "{DB_DIR}": "/run/db", "{FRAGPIPE}": "/fp/fragpipe-24.0"})
        step = root.find("step")
        self.assertEqual((step.get("class"), step.get("name"), step.find("workdir").text), ("CmdMSBooster", "MSBooster", "{WORKDIR}"))
        self.assertEqual([a.text for a in step.find("argv")][-1], "{WORKDIR}/msbooster_params.txt")
        # a file the step's argv names is the step's; FragPipe's run files and the rest are the plan's
        self.assertEqual([(c.get("path"), c.text) for c in step.findall("config")], [("{WORKDIR}/msbooster_params.txt", "useRT = true\n")])
        self.assertEqual(sorted(c.get("path") for c in root.findall("config")),
                         ["{WORKDIR}/fragpipe-files.fp-manifest", "{WORKDIR}/fragpipe.job", "{WORKDIR}/fragpipe.workflow", "{WORKDIR}/g/filelist.txt"])

    def test_the_licence_is_redacted_in_the_argv(self):
        step = ("IonQuant [Work dir: /run/fp_work]", "java -jar /fp/fragpipe-24.0/tools/IonQuant.jar --license /fp/license.dat --threads 4")
        with tempfile.TemporaryDirectory() as tmp:
            root = ptc.load_workdir(workdir(tmp, step))
        self.assertEqual([a.text for a in root.find("step").find("argv")][3:5], ["--license", "{LICENSE}"])

    def test_a_config_naming_the_licence_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ptc.PlanError, "names the licence"):
                ptc.load_workdir(workdir(tmp, MSBOOSTER_STEP, files={"msbooster_params.txt": "license = /fp/license.dat\n"}))

    def argv_of(self, step, **kw):
        with tempfile.TemporaryDirectory() as tmp:
            return [a.text for a in ptc.load_workdir(workdir(tmp, step, **kw)).find("step").find("argv")]

    def test_the_host_name_in_fragpipes_quotes_resolves_from_the_saved_workflow(self):
        # FP-Meta's "Homo sapiens" (CmdMetaproteomics.java:102) prints as ""Homo sapiens"", which the two elements ""Homo and sapiens"" print as too
        step = ("Metaproteomics [Work dir: /run/fp_work]", 'java -cp x.jar DbOptimizer --hostName ""Homo sapiens"" --iterations 1')
        self.assertEqual(self.argv_of(step, workflow="metaproteomics.host-name=Homo sapiens\n")[5], '"Homo sapiens"')

    def test_a_shepherd_config_path_with_a_space_resolves_from_the_work_dir(self):
        # CmdPtmshepherd.java:163 wraps the config path in quotes
        step = ("PTMShepherd [Work dir: /run/fp_work]", 'java -cp x.jar edu.umich.andykong.ptmshepherd.PTMShepherd ""/run/fp_work/a b/shepherd.config""')
        self.assertEqual(self.argv_of(step, files={"a b/shepherd.config": "threads = 4\n"})[-1], '"{WORKDIR}/a b/shepherd.config"')

    def test_cmd_opts_with_quotes_at_word_edges_resolve_from_the_saved_workflow(self):
        # diann.cmd-opts is one argv element (CmdDiann.java:385-387); x" "y in it would end a quoted element early
        step = ("DIA-Quant run DIA-NN [Work dir: /run/fp_work]", 'diann --qvalue 0.01 "--z x" "y" --cfg f')
        self.assertEqual(self.argv_of(step, workflow='diann.cmd-opts=--z x" "y\n')[3], '--z x" "y')

    def test_a_reading_no_run_value_explains_stops(self):
        step = ("Metaproteomics [Work dir: /run/fp_work]", 'java -cp x.jar DbOptimizer --hostName ""Mus musculus"" --iterations 1')
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ptc.PlanError, r"step 1 \(Metaproteomics\): argument 6 reads as '\"\"Mus' or '\"Mus musculus\"', and the run has no value for either"):
                ptc.load_workdir(workdir(tmp, step, workflow="metaproteomics.host-name=Homo sapiens\n"))

    def test_a_line_no_argv_of_the_run_prints_as_stops(self):
        # a spaced argument with a quote at an inner word edge that is no run value: "--a "b" --c" could also end after "b"
        step = ("DIA-Quant run DIA-NN [Work dir: /run/fp_work]", 'diann --qvalue 0.01 "--a "b" --c"')
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ptc.PlanError, r"step 1 \(DIA-Quant run DIA-NN\): no argv FragPipe builds from the run's values prints as"):
                ptc.load_workdir(workdir(tmp, step))

    def test_a_masked_api_key_stops(self):
        step = ("Transfer Learning [Work dir: /run/fp_work]", "python3 -m koina --api-key ******** --x 1")
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ptc.PlanError, r"step 1 \(Transfer Learning\): --api-key is masked"):
                ptc.load_workdir(workdir(tmp, step))

    def test_tools_outside_the_fragpipe_root_get_their_own_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = ptc.load_workdir(workdir(tmp, MSBOOSTER_STEP, tools="/opt/fp-tools"))
        roots = {r.get("name"): r.get("path") for r in root.find("roots")}
        self.assertEqual((roots["{FRAGPIPE}"], roots["{TOOLS}"]), ("/opt", "/opt/fp-tools"))

    def test_a_work_dir_without_one_run_log_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            wd = workdir(tmp, MSBOOSTER_STEP)
            (wd / "log_2026-10-12_11-00-00.txt").write_text("")
            with self.assertRaisesRegex(ptc.PlanError, "expected one FragPipe run log"):
                ptc.load_workdir(wd)


# The goldens of the planner gate: FragPipe work dirs written by FRAGPIPEPLAN (golden*.nf.test) from each
# golden run's saved workflow and manifest; the data types are the manifest's.
GOLDENS = {
    "LFQ-MBR-paper_perfile": "DDA",
    "DIA_SpecLib_Quant_perfile": "DIA",
    "glyco-N-LFQ_perfile": "DDA",
    "DIA_SpecLib_Quant_diaPASEF_perfile": "DIA",
    "DIA_DIA-Umpire_SpecLib_Quant_perfile": "DIA",
    "glyco-O-Pair_perfile": "DDA",
    # MBG
    "glyco-N-LFQ-MBG_perfile": "DDA",
    # the bundled workflows (assets/LFQ-MBR.workflow, assets/TMT6.workflow) on the test_dda_lfq / test_tmtq data
    "LFQ-MBR_test_dda_lfq": "DDA",
    "TMT6_test_tmtq": "DDA",
}
# --threads the pipeline gives the parser (PARSEFRAGPIPEWORKFLOW, params.fragpipe_threads) for workflows without
# workflow.threads; the plans were exported with the same --threads
PARSER_THREADS = {"LFQ-MBR_test_dda_lfq": "2", "TMT6_test_tmtq": "2"}
# Planner vs parser differences, each with its cause; any other difference fails the test.
EXPECTED_DIFFERENCES = {
    # FragPipe writes database_name and num_threads from its run inputs (the plan task's staged FASTA and
    # CPUs); the MSFRAGGER module sets both itself (setFragpipeParamValues, modules/local/msfragger/main.nf:52).
    # In a DIA-only run the parser's args is its DDA base params, which nothing reads; the planner's is the
    # fragger_dia.params FragPipe wrote (configs/ carries the per-data-type files that are read).
    ("msfragger", "args"),
    # FragPipe writes params only for the data types in the manifest; nothing reads the other types
    ("msfragger", "calibrate_mass"),
    # GUI-tab key with no FragPipe step; its effect (IonQuant runs) is in the plan
    ("quantitation", "args"),
    # FragPipe already ANDs useIM / findBestImModel with "the run has timsTOF data" (CmdMSBooster.java:151-157);
    # the parser gives the workflow value and msboosterParamsText ANDs it again; keys are re-ordered there
    ("msbooster", "args"),
}
PARAMS_LINES_SET_BY_MODULE = ("database_name", "num_threads")
MSBOOSTER_ANDED_WITH_TIMSTOF = ("useIM", "findBestImModel")


def strip_module_lines(text):
    return [line for line in text.splitlines() if not line.startswith(PARAMS_LINES_SET_BY_MODULE)]


def msbooster_lines(text, tims_tof=False):
    """MSBooster key = value lines as a set, the timsTOF-gated flags ANDed as msboosterParamsText does."""
    out = set()
    for line in text.splitlines():
        key, _, value = (part.strip() for part in line.partition("="))
        if key in MSBOOSTER_ANDED_WITH_TIMSTOF:
            value = str(value == "true" and tims_tof).lower()
        out.add(f"{key} = {value}")
    return out


def by_class(root, cls):
    return [s for s in root.iter("step") if s.get("class") == cls]


class PlannerParityTest(unittest.TestCase):

    def run_parser(self, name, data_types, out):
        threads = ["--threads", PARSER_THREADS[name]] if name in PARSER_THREADS else []
        subprocess.run([sys.executable, str(PARSER), "--workflow", str(PLANS / f"{name}.workflow"), "--data-types", data_types,
                        "--output-json", "--outdir", out, *threads, *PARSER_TEMPLATES], check=True, capture_output=True)
        subprocess.run([sys.executable, str(PARSER), "--workflow", str(PLANS / f"{name}.workflow"), "--data-types", data_types, *threads,
                        "--outdir", str(Path(out) / "configs"), *PARSER_TEMPLATES], check=True, capture_output=True)
        return json.loads((Path(out) / "tool_configs.json").read_text())

    def test_planner_equals_parser_on_every_golden(self):
        for name, data_types in GOLDENS.items():
            with self.subTest(golden=name), tempfile.TemporaryDirectory() as out:
                planner, files = ptc.tool_configs(ptc.load_workdir(PLANS / name))
                parser = self.run_parser(name, data_types, out)
                self.assertEqual(set(planner), set(parser))
                found = set()
                groups = len(by_class(ptc.load_workdir(PLANS / name), "CmdPhilosopherFilter"))
                for tool in parser:
                    if tool == "abacus" and parser[tool]["run"] and not planner[tool]["run"]:
                        # the parser leaves out FragPipe's "more than one experiment group" condition
                        # (FragpipeRun.java:1782-1809); FRAGPIPE_INFERENCE adds it (ordered.size() > 1)
                        self.assertEqual(groups, 1, f"{name}: abacus off in the plan with {groups} groups")
                        continue
                    self.assertEqual(planner[tool]["run"], parser[tool]["run"], f"{tool}.run")
                    self.assertEqual(planner[tool]["config_type"], parser[tool]["config_type"], f"{tool}.config_type")
                    if not parser[tool]["run"]:
                        continue  # FragPipe builds no step for a disabled tool, and nothing reads its args
                    found |= {(tool, f) for f in set(planner[tool]) | set(parser[tool]) if planner[tool].get(f) != parser[tool].get(f)}
                self.assertLessEqual(found, EXPECTED_DIFFERENCES, "unexplained planner/parser differences")
                if ("msfragger", "calibrate_mass") in found:
                    a, b = planner["msfragger"]["calibrate_mass"], parser["msfragger"]["calibrate_mass"]
                    self.assertEqual(a, {k: v for k, v in b.items() if k in a})
                if ("msbooster", "args") in found:
                    self.assertEqual(msbooster_lines(planner["msbooster"]["args"]), msbooster_lines(parser["msbooster"]["args"]))
                # the legacy configs/ files the subworkflows read: per-data-type MSFragger params, MSBooster params
                for file_name, text in files.items():
                    if file_name == "msbooster.config" and not parser["msbooster"]["run"]:
                        self.assertEqual(text, "", "MSBooster does not run: empty msbooster.config")
                        continue
                    parser_text = (Path(out) / "configs" / file_name).read_text()
                    if file_name.endswith(".params"):
                        self.assertEqual(strip_module_lines(text), strip_module_lines(parser_text), file_name)
                    else:
                        self.assertEqual(msbooster_lines(text), msbooster_lines(parser_text), file_name)


PHI = "{FRAGPIPE}/tools/Philosopher/philosopher-v5.1.3-RC9"


# The committed work dirs name FragPipe's install as @ROOT@/fragpipe-24.0 with the academic jar names; a real
# dry run names the image's install path, and another FragPipe edition names its jars <Tool>-<Edition>-<version>.jar
INSTALL = "/opt/fragpipe/fragpipe-24.0"
EDITION_JARS = re.compile(r"\b(MSFragger|IonQuant|diaTracer)-(?=[0-9])")


def as_installed(src, dst, edition):
    """A copy of a fixture work dir with FragPipe at INSTALL, and edition jar names (<Tool>-Edition-<version>.jar) when ``edition``."""
    shutil.copytree(src, dst)
    for path in dst.rglob("*"):
        if path.is_file():
            text = path.read_text(encoding="utf-8")
            text = text.replace("@ROOT@/fragpipe-24.0", INSTALL)
            if edition:
                text = EDITION_JARS.sub(lambda m: f"{m.group(1)}-Edition-", text)
            path.write_text(text, encoding="utf-8")


class InstallLayoutTest(unittest.TestCase):
    """The tool configs do not depend on FragPipe's install path or on the jar names with or without an edition word."""

    def test_install_path_and_jar_names_give_the_same_tool_configs(self):
        for name in GOLDENS:
            expected = ptc.tool_configs(ptc.load_workdir(PLANS / name))
            for edition in (False, True):
                with self.subTest(golden=name, edition=edition), tempfile.TemporaryDirectory() as tmp:
                    wd = Path(tmp) / name
                    as_installed(PLANS / name, wd, edition)
                    log = next(wd.glob("log_*.txt")).read_text(encoding="utf-8")
                    self.assertIn(INSTALL, log)
                    self.assertEqual("-Edition-" in log, edition)
                    self.assertEqual(ptc.tool_configs(ptc.load_workdir(wd)), expected)


class QuantAndProphetExtractorsTest(unittest.TestCase):
    """FreeQuant, LabelQuant, PTMProphet, iProphet, MBG and the Abacus TMT label (plan argv as FragPipe builds it)."""

    def test_tmt_freequant_does_not_run_label_free_freequant(self):
        fq = step(1, "CmdFreequant", "FreeQuant", [PHI, "freequant", "--ptw", "0.4", "--tol", "10", "--dir", "/w/lcms"])
        lq = step(2, "CmdLabelquant", "Quant (Isobaric)", [PHI, "labelquant", "--tol", "10", "--plex", "18", "--annot", "/w/P1_annotation.txt",
                                                           "--brand", "tmt", "--raw", "--dir", "/w/lcms"])
        configs, _files = ptc.tool_configs(plan(fq, lq))
        self.assertEqual(configs["freequant"], {"run": False, "args": "--ptw 0.4 --tol 10", "config_type": "cli"})
        self.assertEqual(configs["labelquant"], {"run": True, "args": "--tol 10 --plex 18", "config_type": "cli", "brand": "tmt"})
        configs, _files = ptc.tool_configs(plan(fq))
        self.assertTrue(configs["freequant"]["run"], "without LabelQuant the FreeQuant step is the label-free one")
        with self.assertRaisesRegex(ptc.PlanError, "2 FreeQuant steps for 1 LabelQuant"):
            ptc.tool_configs(plan(fq, step(3, "CmdFreequant", "FreeQuant", [PHI, "freequant", "--ptw", "0.4", "--tol", "10", "--dir", "/w/lcms"]), lq))

    def test_labelquant_groups_with_different_brands_stop(self):
        def lq(o, brand):
            return step(o, "CmdLabelquant", "Quant (Isobaric)", [PHI, "labelquant", "--plex", "6", "--annot", "a.txt", "--brand", brand, "--dir", "/w"])
        with self.assertRaisesRegex(ptc.PlanError, "different --brand"):
            ptc.tool_configs(plan(lq(1, "tmt"), lq(2, "itraq")))

    def test_ptmprophet_drops_fragpipes_maxthreads_and_files(self):
        pp = step(1, "CmdPtmProphet", "PTMProphet", ["/f/PTMProphetParser", "NOSTACK", "STY:79.96633", "MINPROB=0.5", "MAXTHREADS=1",
                                                     "interact-a.pep.xml", "interact-a.mod.pep.xml"])
        configs, _files = ptc.tool_configs(plan(pp))
        self.assertEqual(configs["ptmprophet"]["args"], "NOSTACK STY:79.96633 MINPROB=0.5")
        with self.assertRaisesRegex(ptc.PlanError, "MAXTHREADS=1"):
            ptc.tool_configs(plan(step(1, "CmdPtmProphet", "PTMProphet", ["/f/PTMProphetParser", "NOSTACK", "a.pep.xml", "a.mod.pep.xml"])))

    def test_iprophet_args_without_threads_and_pepxmls(self):
        ip = step(1, "CmdIprophet", "iProphet", [PHI, "iprophet", "--decoy", "rev_", "--nonsp", "--output", "combined", "--threads", "24",
                                                 "/w/A01/interact-a.pep.xml", "/w/B01/interact-b.pep.xml"])
        configs, _files = ptc.tool_configs(plan(ip))
        self.assertEqual(configs["iprophet"]["args"], "--decoy rev_ --nonsp --output combined")

    def test_mbg_runtime_values_become_tokens_and_im_ms_comes_from_the_workflow(self):
        def mbg(o, g):
            return step(o, "CmdMBGMatch", "MBG", ["java", "-cp", "m.jar", "com.mbg.MBG", "--match", "--psm", f"/w/{g}/psm.tsv",
                                                  "--manifest", "/w/fragpipe-files.fp-manifest", "--residuedb",
                                                  "{FRAGPIPE}/lib/../tools/Glycan_Databases/glycan_residues.txt",
                                                  "--nopasef", "true", "--numthreads", "48", "--runtmt", "false"])
        configs, _files = ptc.tool_configs(plan(mbg(1, "cancer"), mbg(2, "saliva"), workflow="workflow.input.data-type.im-ms=true\n"))
        self.assertEqual(configs["mbg"]["args"], ["--match", "--psm", "@PSM@", "--manifest", "@MANIFEST@", "--residuedb",
                                                  "@TOOLS_DIR@/Glycan_Databases/glycan_residues.txt", "--nopasef", "@NOPASEF@",
                                                  "--numthreads", "@THREADS@", "--runtmt", "false"])
        self.assertTrue(configs["mbg"]["im_ms"])

    def test_abacus_tmt_label_from_labels_plex_and_the_workflow(self):
        ab = step(1, "CmdPhilosopherAbacus", "PhilosopherAbacus", [PHI, "abacus", "--razor", "--reprint", "--tag", "rev_", "--protein",
                                                                   "--labels", "--plex", "16", "G1", "G2"])
        configs, _files = ptc.tool_configs(plan(ab, workflow="tmtintegrator.channel_num=TMT-16\n"))
        self.assertEqual(configs["abacus"], {"run": True, "args": "--razor --reprint --tag rev_ --protein", "config_type": "cli",
                                             "labels_channel_num": "TMT-16"})


if __name__ == "__main__":
    unittest.main()
