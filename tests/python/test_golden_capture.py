"""Unit tests for tests/tools/golden_capture.py and tests/tools/golden_step_check.py."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import fragpipe_log_to_manifest as fpm  # noqa: E402
import golden_capture as cap  # noqa: E402
import golden_step_check as gsc  # noqa: E402

# golden_step_check placeholder of the FragPipe install (a plain string: nf-core lint reads escaped f-string braces as Jinja)
FRAGPIPE_TOKEN = "{FRAGPIPE}"



class TestLicenseScan(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.f = Path(self._tmp.name) / "x.txt"

    def tearDown(self):
        self._tmp.cleanup()

    def test_variants_detected(self):
        for text in ("a --license b", "a -Licence b", "LICENSE_FILE=x", "licence-key", "/a/Licence.DAT"):
            with self.subTest(text=text):
                self.f.write_text(text)
                self.assertTrue(cap._mentions_license(self.f))

    def test_binary_detection(self):
        d = Path(self._tmp.name)
        raw = d / "a.raw"
        raw.write_bytes(b"license.dat --license")
        self.assertTrue(cap._is_binary(raw))
        nul = d / "x.dat"
        nul.write_bytes(b"abc\x00license.dat")
        self.assertTrue(cap._is_binary(nul))
        (d / "run.d").mkdir()
        tdf = d / "run.d" / "analysis.tdf"
        tdf.write_bytes(b"x")
        self.assertTrue(cap._is_binary(tdf))
        txt = d / "fragpipe.job"
        txt.write_text("--license x")
        self.assertFalse(cap._is_binary(txt))
        mz = d / "a.mzML"
        mz.write_text("<mzML/>")
        self.assertFalse(cap._is_binary(mz))

    def test_protein_text_not_detected(self):
        self.f.write_text("DNA replication licensing factor MCM2\n")
        self.assertFalse(cap._mentions_license(self.f))

    def test_match_in_a_long_line_is_found(self):
        self.f.write_bytes(b"x" * (5 << 20) + b" license-path=/a\n")
        self.assertTrue(cap._mentions_license(self.f))

    def test_redaction_placeholders_are_accepted(self):
        for text in (
            "java -jar MSFragger.jar p --license {LICENSE} a.mzML\n",
            "ionquant --license={LICENSE} --threads 4\n",
            "java -Dk={LICENSE_PATH} -jar x\n",
            "Using license {LICENSE_PATH}\n".replace("license ", "--license "),
        ):
            with self.subTest(text=text):
                self.f.write_text(text)
                self.assertFalse(cap._mentions_license(self.f))

    def test_real_values_next_to_placeholders_still_rejected(self):
        for text in (
            "x --license {LICENSE} --license /opt/key.lic\n",
            "x {LICENSE_PATH} /opt/fp/license.dat\n",
            "x --license=abc\n",
            "x --license {LICENSEX}\n",
            "Licence file: {LICENSE}\n",
        ):
            with self.subTest(text=text):
                self.f.write_text(text)
                self.assertTrue(cap._mentions_license(self.f))

class TestArgvNormalisation(unittest.TestCase):
    """golden_step_check argv normalisation, without a captured run."""

    def test_filemove_commons_io_classpath(self):
        # ToolingUtils.java:144-148: FileMove's second -cp entry is the jar commons-io FileUtils was loaded from
        fp = "{FRAGPIPE}/lib/fragpipe-24.0.jar"
        norm = ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/{gsc.COMMONS_IO_TOKEN}", "org.nesvilab.utils.FileMove", "--no-err", "a", "b"]
        for jar in ("commons-io-2.18.0.jar", "batmass-io-1.36.5.jar", "batmass-io-1.36.6.jar"):
            argv = ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/{jar}", "org.nesvilab.utils.FileMove", "--no-err", "a", "b"]
            self.assertEqual(gsc.accepted_normalisations(argv), norm, jar)
        exp = ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/batmass-io-1.36.5.jar", "org.nesvilab.utils.FileMove", "a", "b"]
        act = ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/commons-io-2.18.0.jar", "org.nesvilab.utils.FileMove", "a", "b"]
        self.assertEqual(gsc.check_argv(exp, act)["status"], "MATCH")

    def test_filemove_classpath_negative(self):
        fp = "{FRAGPIPE}/lib/fragpipe-24.0.jar"
        # another jar in FileMove's -cp, the same jars for another class, and the destination or a third entry stay as they are
        for argv in (
            ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/guava-33.jar", "org.nesvilab.utils.FileMove", "a", "b"],
            ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/batmass-io-1.36.6.jar", "org.nesvilab.utils.FileMoveWithExtension", "a"],
            ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/commons-io-2.18.0.jar", "org.nesvilab.fragpipe.tools.diann.Propagation", "a"],
            ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/commons-io-2.18.0.jar:{FRAGPIPE_TOKEN}/lib/x.jar", "org.nesvilab.utils.FileMove", "a"],
            ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/commons-io-2.17.0.jar", "org.nesvilab.utils.FileMove", "a", "b"],
        ):
            self.assertEqual(gsc.accepted_normalisations(argv), argv, argv)
        exp = ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/batmass-io-1.36.5.jar", "org.nesvilab.utils.FileMove", "a", "b"]
        act = ["java", "-cp", f"{fp}:{FRAGPIPE_TOKEN}/lib/commons-io-2.18.0.jar", "org.nesvilab.utils.FileMove", "a", "c"]
        self.assertEqual(gsc.check_argv(exp, act)["status"], "DIFF")

    def test_accepted_normalisations(self):
        argv = ["java", "-Duser.home=/t/x", "-Djava.io.tmpdir=/t/y", "-Xmx62G", "-jar", "IonQuant.jar", "--license", "/opt/fp/license.dat", "--threads", "4"]
        self.assertEqual(
            gsc.accepted_normalisations(argv),
            ["java", "-Xmx<N>", "-jar", "IonQuant.jar", "--license", "{LICENSE}", "--threads", "4"],
        )
        self.assertEqual(gsc.accepted_normalisations(["--license=/x/license.dat"]), ["--license={LICENSE}"])

    def test_other_java_options_are_not_normalised(self):
        self.assertEqual(gsc.accepted_normalisations(["-Dfile.encoding=UTF-8", "-Xms4G"]), ["-Dfile.encoding=UTF-8", "-Xms4G"])

    def test_actual_lcms_dir_and_workdir_map_to_golden_tokens(self):
        step = {"cwd": "{WORKDIR}/P1", "inputs": [], "outputs": []}
        argv = gsc.actual_argv("tool --dir /t/task/lcms_files /t/task/g.tsv", step, "/t/task/P1", "/opt/fp", "/t/task", "/t/task/lcms_files")
        self.assertEqual(argv, ["tool", "--dir", "{LCMS_DIR}", "{WORKDIR}/g.tsv"])

    def test_actual_db_dir_is_a_module_path_prefix(self):
        args = SimpleNamespace(actual_cwd="/t/task", actual_workdir=None, actual_db_dir="/g/db")
        self.assertEqual(gsc.actual_prefixes(args, ["/fp/work"]), ["/t/task", "/fp/work", "/g/db"])
        side = gsc.gc.Side(prefixes=gsc.actual_prefixes(args, []))
        self.assertEqual(side.rewrite_paths(b'local_path="/g/db/x.fasta"'), b'local_path="' + gsc.gc.RUN_TOKEN + b'/x.fasta"')

    def test_without_actual_db_dir_prefixes_are_unchanged(self):
        args = SimpleNamespace(actual_cwd="/t/task", actual_workdir="/t/wd")
        self.assertEqual(gsc.actual_prefixes(args, ["/fp/work"]), ["/t/task", "/fp/work", "/t/wd"])

    def test_staged_run_lcms_file_takes_its_input_token_when_the_step_records_no_inputs(self):
        step = {"cwd": None, "inputs": [], "outputs": [], "argv": []}
        run_inputs = [{"path": "{INPUT:run.mzML}", "experiment": "", "data_type": "DIA"}]
        got = gsc.actual_argv("java -cp x.jar Main run.mzML /task/umpire-se.params", step, "/nowhere", "/opt/fp", "/task", None, run_inputs)
        self.assertEqual(got, ["java", "-cp", "x.jar", "Main", "{INPUT:run.mzML}", "{WORKDIR}/umpire-se.params"])

    def test_staged_database_in_the_input_dir_takes_its_input_token(self):
        step = {"cwd": "{WORKDIR}", "inputs": [], "outputs": [], "argv": []}
        run_inputs = [{"path": "{INPUT:db.fasta}"}]
        got = gsc.actual_argv("java -jar m.jar --projectDir /task/project --fastaFile /task/db.fasta --outdir /task/project", step, "/task/project", "/opt/fp", "/task/project", None, run_inputs, "/task")
        self.assertEqual(got, ["java", "-jar", "m.jar", "--projectDir", "{CWD}", "--fastaFile", "{INPUT:db.fasta}", "--outdir", "{CWD}"])

    def test_module_argv_split_keeps_empty_and_quoted_args(self):
        step = {"cwd": "{WORKDIR}/C02", "inputs": [], "outputs": [], "argv": []}
        got = gsc.actual_argv('/opt/fp/tools/x --mod-tag "" --name "a b" -v', step, "/task", "/opt/fp")
        self.assertEqual(got, ["{FRAGPIPE}/tools/x", "--mod-tag", "", "--name", "a b", "-v"])

    def test_fragpipe_root_from_tool_path(self):
        # executable, -jar, -cp (classpath entries) and FragPipe's lib/* all name the install root
        self.assertEqual(gsc.fragpipe_root("/fp/fragpipe-24/tools/Philosopher/philosopher-v5 filter --pepxml /task/C02"), "/fp/fragpipe-24")
        self.assertEqual(gsc.fragpipe_root("java -Xmx6G -jar /fp/fragpipe-24/tools/TMT-Integrator-6.1.3.jar conf.yml /task/psm.tsv"), "/fp/fragpipe-24")
        self.assertEqual(gsc.fragpipe_root('java -cp "/fp/fragpipe-24/lib/*" org.nesvilab.fragpipe.tools.x /task/a'), "/fp/fragpipe-24")
        self.assertEqual(gsc.fragpipe_root("java -cp /x/a.jar:/fp/fragpipe-24/tools/batmass-io.jar Main /task/a"), "/fp/fragpipe-24")

    def test_fragpipe_root_skips_other_lib_dirs(self):
        # an absolute java under /usr/lib/jvm is not the FragPipe root
        cmd = "/usr/lib/jvm/java-17-openjdk-amd64/bin/java -cp /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/x.jar Main /task/a"
        self.assertEqual(gsc.fragpipe_root(cmd), "/fragpipe_bin/fragpipe-24.0/fragpipe-24.0")
        self.assertEqual(gsc.fragpipe_root("/opt/tools/fragpipe-24.0/tools/MSFragger/msfragger -v"), "/opt/tools/fragpipe-24.0")

    def test_fragpipe_root_two_roots_fail(self):
        with self.assertRaises(ValueError):
            gsc.fragpipe_root("java -cp /a/fragpipe-23/lib/*:/b/fragpipe-24/tools/x.jar Main")

    def test_fragpipe_root_without_tool_path_fails(self):
        with self.assertRaises(ValueError):
            gsc.fragpipe_root("tool --dir /t/task/lcms_files")

    def test_actual_argv_derives_fragpipe_root(self):
        step = {"inputs": [], "outputs": [], "cwd": "{WORKDIR}"}
        self.assertEqual(gsc.actual_argv("/fp/fragpipe-24/tools/x -v", step, "/task", None)[0], "{FRAGPIPE}/tools/x")

    def test_literal_quotes_in_an_argument_are_compared_exactly(self):
        # CmdPtmshepherd.java:163 passes '"<config>"' (literal quotes) through a ProcessBuilder; FragPipe's log
        # prints it as is (only whitespace arguments get an extra wrap), the module log is printf %q.
        # The golden side goes through split_command, so this fails when split_command strips the quotes.
        step = {"cwd": "{WORKDIR}", "inputs": [], "outputs": [], "argv": []}
        golden = fpm.Normalizer({"/w": "{WORKDIR}"}).argv(fpm.split_command('java Main "/w/shepherd.config"'))
        quoted = gsc.actual_argv(r"java Main \"/w/shepherd.config\"", step, "/nowhere", "/opt/fp", "/w")
        bare = gsc.actual_argv("java Main /w/shepherd.config", step, "/nowhere", "/opt/fp", "/w")
        self.assertEqual(gsc.check_argv(golden, quoted)["status"], "MATCH")
        self.assertEqual(gsc.check_argv(golden, bare)["status"], "DIFF")



class TestWorkspaceSteps(unittest.TestCase):
    def steps(self, next_status):
        return [
            {"step_index": 23, "cmd_name": "PhilosopherDbAnnotate", "cwd": "{WORKDIR}/C02", "argv": {"status": "MATCH"}, "files": []},
            {"step_index": 24, "cmd_name": "PhilosopherFilter", "cwd": "{WORKDIR}/C02", "argv": {"status": "MATCH"}, "files": [{"status": next_status}]},
            {"step_index": 25, "cmd_name": "PhilosopherDbAnnotate", "cwd": "{WORKDIR}/D02", "argv": {"status": "MATCH"}, "files": []},
        ]

    def test_workspace_step_is_verified_by_the_matching_next_step_in_its_dir(self):
        st = self.steps("MATCH")
        gsc.resolve_workspace_steps(st)
        self.assertEqual([f["status"] for f in st[0]["files"]], ["VIA_CONSUMER"])
        self.assertEqual(st[2]["files"], [], "no later step in D02: stays unverified")

    def test_differing_next_step_leaves_the_workspace_step_unverified(self):
        st = self.steps("DIFF")
        gsc.resolve_workspace_steps(st)
        self.assertEqual(st[0]["files"], [])



class TestEditionJars(unittest.TestCase):
    def test_edition_word_of_licensed_jars_is_not_compared(self):
        # FragPipe finds MSFragger/IonQuant/diaTracer under either name (TabConfig.java:176-178, :188-209)
        argv = ["java", "-jar", "{FRAGPIPE}/tools/MSFragger-4.4.1/MSFragger-Plus-4.4.1.jar", "-cp", "/x/IonQuant-Plus-1.11.20.jar:/x/diaTracer-Plus-2.2.1.jar"]
        self.assertEqual(
            gsc.accepted_normalisations(argv),
            ["java", "-jar", "{FRAGPIPE}/tools/MSFragger-4.4.1/MSFragger-4.4.1.jar", "-cp", "/x/IonQuant-1.11.20.jar:/x/diaTracer-2.2.1.jar"],
        )
        self.assertEqual(gsc.check_argv(["java", "-jar", "/f/MSFragger-4.4.1.jar"], ["java", "-jar", "/f/MSFragger-Plus-4.4.1.jar"])["status"], "MATCH")

    def test_other_jars_and_versions_are_compared(self):
        self.assertEqual(gsc.accepted_normalisations(["/x/MSBooster-1.4.17.jar", "/x/MSFragger-4.4.2.jar"]), ["/x/MSBooster-1.4.17.jar", "/x/MSFragger-4.4.2.jar"])
        self.assertEqual(gsc.check_argv(["java", "-jar", "/f/MSFragger-4.4.1.jar"], ["java", "-jar", "/f/MSFragger-Plus-4.4.2.jar"])["status"], "DIFF")


if __name__ == "__main__":
    unittest.main()
