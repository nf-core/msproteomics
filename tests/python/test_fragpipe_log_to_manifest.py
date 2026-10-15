"""Unit tests for tests/tools/fragpipe_log_to_manifest.py (run: make test-python)."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import fragpipe_log_to_manifest as fpm  # noqa: E402


class TestRedaction(unittest.TestCase):
    def test_license_argument_value_is_redacted(self):
        line = "java -jar /x/MSFragger.jar /w/fragger.params --license /secret/path/lic.dat /in/a.mzML"
        self.assertEqual(
            fpm.redact_line(line), "java -jar /x/MSFragger.jar /w/fragger.params --license {LICENSE} /in/a.mzML"
        )

    def test_license_path_token_replaced_command_kept(self):
        line = "java -Dx=1 -cp /t/IonQuant.jar ionquant.IonQuant --threads 4 /opt/fp/license.dat --mbr 0"
        self.assertEqual(
            fpm.redact_line(line), "java -Dx=1 -cp /t/IonQuant.jar ionquant.IonQuant --threads 4 {LICENSE_PATH} --mbr 0"
        )

    def test_already_scrubbed_line_is_unchanged(self):
        line = "java -jar /x/MSFragger.jar p --license {LICENSE} {LICENSE_PATH} a.mzML"
        self.assertEqual(fpm.redact_line(line), line)
        self.assertFalse(fpm.contains_license_material(line))

    def test_license_token_in_command_is_redacted(self):
        self.assertEqual(fpm.redact_line("/x/tool --key=/a/license.dat -v"), "/x/tool --key={LICENSE_PATH} -v")
        self.assertEqual(fpm.redact_line("/x/tool --licensekey=abc -v"), "/x/tool --licensekey={LICENSE} -v")

    def test_license_stdout_line_is_dropped(self):
        self.assertIsNone(fpm.redact_line("Licensed to Some Org until 2027"))
        self.assertNotIn("Licensed", fpm.redact_text("a\nLicensed to X\nb"))

    def test_lines_without_license_are_unchanged(self):
        line = "Process 'MSFragger' finished, exit code: 0"
        self.assertEqual(fpm.redact_line(line), line)


class TestRedactionVariants(unittest.TestCase):
    """Every flag spelling/separator; argv shape kept; no command line dropped."""

    CASES = {
        "java -jar M.jar p --license /s/k.dat a.mzML": "java -jar M.jar p --license {LICENSE} a.mzML",
        "/t/tool --License KEY y": "/t/tool --License {LICENSE} y",
        "/t/tool --LICENSE KEY y": "/t/tool --LICENSE {LICENSE} y",
        "/t/tool --licence KEY y": "/t/tool --licence {LICENSE} y",
        "/t/tool -license KEY y": "/t/tool -license {LICENSE} y",
        "/t/tool ---licence KEY y": "/t/tool ---licence {LICENSE} y",
        "/t/tool --license=KEY y": "/t/tool --license={LICENSE} y",
        "/t/tool --Licence=KEY y": "/t/tool --Licence={LICENSE} y",
        "/t/tool --license= y": "/t/tool --license={LICENSE} y",
        "/t/tool --license-key KEY y": "/t/tool --license-key {LICENSE} y",
        "/t/tool --licensefile=/a/b y": "/t/tool --licensefile={LICENSE} y",
        "/t/tool -Dk=/opt/fp/license.dat -jar x": "/t/tool -Dk={LICENSE_PATH} -jar x",
        "/t/tool /opt/fp/LICENCE.DAT y": "/t/tool {LICENSE_PATH} y",
        "/t/tool --license --threads 4": "/t/tool --license {LICENSE} 4",
        "/t/tool --license {LICENSE} y": "/t/tool --license {LICENSE} y",
        "/t/tool --license {LICENSE_PATH} y": "/t/tool --license {LICENSE_PATH} y",
    }

    def test_flag_variants(self):
        for line, want in self.CASES.items():
            with self.subTest(line=line):
                self.assertEqual(fpm.redact_line(line), want)
                self.assertEqual(len(fpm.redact_line(line).split(" ")), len(line.split(" ")), "argv shape changed")

    def test_empty_value_keeps_shape_and_redacts_next(self):
        line = "/t/tool --license  KEY --x 1"
        out = fpm.redact_line(line)
        self.assertEqual(out, "/t/tool --license  {LICENSE} --x 1")
        self.assertEqual(fpm.split_command(out), ["/t/tool", "--license", "", "{LICENSE}", "--x", "1"])

    def test_value_starting_with_a_dash_is_redacted(self):
        self.assertEqual(fpm.redact_line("/t/tool --license -K3Y/abc --mbr 0"), "/t/tool --license {LICENSE} --mbr 0")
        self.assertEqual(fpm.redact_line("/t/tool --license --mbr 0"), "/t/tool --license {LICENSE} 0")

    def test_commands_with_any_start_are_kept(self):
        for line in (
            "python3 s.py --license-key abc",
            "bash run.sh --license abc",
            "./bin/tool --license abc",
            "../tool --license abc",
            "bin/tool --license abc",
            "tool.py --license abc",
            "diann-1.8 --license abc --threads 2",
            "mytool --licence abc --threads 2",
        ):
            with self.subTest(line=line):
                out = fpm.redact_line(line)
                self.assertIsNotNone(out)
                self.assertNotIn("abc", out)
                self.assertEqual(len(out.split(" ")), len(line.split(" ")))

    def test_line_after_step_header_is_a_command_whatever_it_looks_like(self):
        out = fpm.redact_lines(["MSFragger [Work dir: /w]", "weird start --license K tail"])
        self.assertEqual(out, ["MSFragger [Work dir: /w]", "weird start --license {LICENSE} tail"])

    def test_informational_lines_are_dropped(self):
        for line in ("Licensed to Some Org until 2027", "Usage: --license /p/key", "LICENCE valid - ok", "[license /p/k redacted]"):
            with self.subTest(line=line):
                self.assertIsNone(fpm.redact_line(line))

    def test_detector_ignores_only_placeholders_and_flag_names(self):
        self.assertFalse(fpm.contains_license_material('"--license", "{LICENSE}", "-licence-file={LICENSE_PATH}"'))
        self.assertTrue(fpm.contains_license_material("Licensed to X"))
        self.assertTrue(fpm.contains_license_material("/opt/license.dat"))


class TestStepKinds(unittest.TestCase):
    def test_all_24_0_cmd_names_known(self):
        for name in ("PairScans", "OPair", "MBG", "FPOP", "FPOP-adapter", "iProphet", "Skyline", "WriteSubMzml",
                     "Quant (Isobaric)", "Transfer Learning", "ExportMatchedFragments", "MSFragger Digest"):
            with self.subTest(name=name):
                self.assertIsNotNone(fpm.step_kind(name))

    def test_sub_step_by_known_prefix(self):
        self.assertEqual(fpm.step_kind("IonQuant: write report"), ("CmdIonquant", "IonQuant: write report"))
        self.assertEqual(fpm.step_kind("DIA-Quant run DIA-NN")[0], "CmdDiann")
        self.assertIsNone(fpm.step_kind("IonQuantX"))
        self.assertIsNone(fpm.step_kind("Totally new"))


class TestSplitCommand(unittest.TestCase):
    def test_plain_tokens(self):
        self.assertEqual(fpm.split_command("a -b c"), ["a", "-b", "c"])

    def test_quoted_argument_with_spaces(self):
        self.assertEqual(fpm.split_command('java -cp "/a b/x.jar" Main'), ["java", "-cp", "/a b/x.jar", "Main"])

    def test_quoted_token_without_space_keeps_its_literal_quotes(self):
        # CmdPtmshepherd.java:163 adds the quotes itself; FragPipe's log wraps only whitespace arguments
        self.assertEqual(fpm.split_command('java Main "/w/shepherd.config"'), ["java", "Main", '"/w/shepherd.config"'])

    def test_wrapped_argument_with_literal_quotes(self):
        self.assertEqual(fpm.split_command('a ""/w x/s.config"" b'), ["a", '"/w x/s.config"', "b"])

    def test_empty_argument_is_preserved(self):
        self.assertEqual(fpm.split_command("fragreporter --mod-tag  --x 1"), ["fragreporter", "--mod-tag", "", "--x", "1"])

    def test_trailing_empty_argument_is_preserved(self):
        self.assertEqual(fpm.split_command("a --mod-tag "), ["a", "--mod-tag", ""])

    def test_double_space_is_an_empty_argument(self):
        self.assertEqual(fpm.split_command("a --mod-tag  --x 1"), ["a", "--mod-tag", "", "--x", "1"])

    def test_unterminated_quote_raises(self):
        with self.assertRaises(fpm.LogFormatError):
            fpm.split_command('a "b c')



class TestMsfraggerOnUmpirePseudoSpectra(unittest.TestCase):
    """DIA-Umpire writes <name>_Q1.mzML next to the input; MSFragger searches them with the DDA params,
    and they are not in the manifest's LC-MS list (golden set DIA_DIA-Umpire_SpecLib_Quant_perfile)."""

    def test_pseudo_spectra_in_argv_are_inputs_and_name_the_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            work.mkdir()
            fasta = Path(d) / "db.fasta"
            fasta.write_text(">sp|P1|A_HUMAN\nPEPTIDEK\n")
            params = work / "fragger_dda.params"
            params.write_text(f"database_name = {fasta}\t\t# Path\ndata_type = 0\noutput_format = pepXML_pin\n")
            q1 = work / "run_Q1.mzML"
            q1.write_text("<mzML/>\n")
            for ext in (".pepXML", ".pin"):
                (work / f"run_Q1{ext}").write_text("x\n")
            step = {
                "step_index": 37, "tool": "MSFragger", "cwd_raw": str(work),
                "argv_raw": ["java", "-jar", "-Dfile.encoding=UTF-8", "-Xmx100G", "/fp/tools/MSFragger.jar", str(params), "--license", "{LICENSE}", str(q1)],
            }
            files = {str(p): {} for p in work.iterdir()}
            files[str(fasta)] = {}
            fpm.derive_io(step, files, {"lcms_files": [{"path": str(Path(d) / "input" / "run.mzML")}]}, {}, {})
            self.assertIn((str(q1), "argv"), step["_inputs"])
            self.assertIn((str(work / "run_Q1.pepXML"), "convention:output_format"), step["_outputs"])
            self.assertIn((str(work / "run_Q1.pin"), "convention:output_format"), step["_outputs"])

    def test_spectra_extensions_match_in_any_case(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            work.mkdir()
            fasta = Path(d) / "db.fasta"
            fasta.write_text(">sp|P1|A_HUMAN\nPEPTIDEK\n")
            params = work / "fragger.params"
            params.write_text(f"database_name = {fasta}\t\t# Path\ndata_type = 0\noutput_format = pepXML\n")
            raw = Path(d) / "RUN.RAW"
            raw.write_text("x\n")
            step = {
                "step_index": 5, "tool": "MSFragger", "cwd_raw": str(work),
                "argv_raw": ["java", "-jar", "/fp/tools/MSFragger.jar", str(params), "--license", "{LICENSE}", str(raw)],
            }
            files = {str(p): {} for p in work.iterdir()}
            files[str(fasta)] = {}
            fpm.derive_io(step, files, {"lcms_files": [{"path": str(raw)}]}, {}, {})
            self.assertIn((str(raw), "argv"), step["_inputs"])

    def test_diatracer_pseudo_spectra_are_msfragger_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            work.mkdir()
            fasta = Path(d) / "db.fasta"
            fasta.write_text(">sp|P1|A_HUMAN\nPEPTIDEK\n")
            params = work / "fragger.params"
            params.write_text(f"database_name = {fasta}\t\t# Path\ndata_type = 0\noutput_format = pepXML\n")
            pseudo = work / "run_diatracer.mzML"
            pseudo.write_text("<mzML/>\n")
            step = {
                "step_index": 6, "tool": "MSFragger", "cwd_raw": str(work),
                "argv_raw": ["java", "-jar", "/fp/tools/MSFragger.jar", str(params), "--license", "{LICENSE}", str(pseudo)],
            }
            files = {str(p): {} for p in work.iterdir()}
            files[str(fasta)] = {}
            fpm.derive_io(step, files, {"lcms_files": [{"path": str(Path(d) / "input" / "run.d")}]}, {}, {})
            self.assertIn((str(pseudo), "argv"), step["_inputs"])

    def test_msfragger_without_any_spectra_in_argv_still_raises(self):
        step = {"step_index": 1, "tool": "MSFragger", "cwd_raw": "/w", "argv_raw": ["java", "-jar", "/fp/MSFragger.jar", "/nonexistent.params"]}
        with self.assertRaises((fpm.LogFormatError, OSError)):
            fpm.derive_io(step, {}, {"lcms_files": []}, {}, {})



GLYCO_LOG = Path(__file__).resolve().parents[1] / "fixtures" / "golden" / "fragpipe24_glyco_lfq_log.trimmed.txt"


class TestParallelGlycoLog(unittest.TestCase):
    """Real FragPipe 24.0 glyco-N-LFQ log: PeptideProphet runs per group in parallel."""

    @classmethod
    def setUpClass(cls):
        cls.text = GLYCO_LOG.read_text(encoding="utf-8")
        cls.parsed = fpm.parse_log(cls.text)
        cls.by = {s["step_index"]: s for s in cls.parsed["steps"]}

    def test_all_29_steps_parsed(self):
        self.assertEqual(len(self.parsed["steps"]), 29)
        self.assertTrue(all(s["exit_code"] == 0 for s in self.parsed["steps"]))

    def test_parallel_peptideprophet_steps_recorded(self):
        pp = [s for s in self.parsed["steps"] if s["name"] == "PeptideProphet"]
        self.assertEqual(len(pp), 2)
        self.assertEqual(pp[0]["concurrent_with"], [pp[1]["step_index"]])
        self.assertEqual(pp[0]["parallel_batch"], pp[1]["parallel_batch"])
        self.assertEqual(pp[0]["finish_attribution_ambiguous_with"], [pp[1]["step_index"]])

    def test_new_step_kinds_recognized(self):
        names = {s["name"] for s in self.parsed["steps"]}
        self.assertTrue({"Rewrite pepxml", "AppendFile", "PTMShepherd", "MSFragger move tsv"} <= names)

    def test_msfragger_command_lines_kept_with_license_placeholder(self):
        mf = [s for s in self.parsed["steps"] if s["name"] == "MSFragger"]
        self.assertEqual(len(mf), 2)
        for s in mf:
            self.assertIn("--license", s["argv_raw"])
            self.assertEqual(s["argv_raw"][s["argv_raw"].index("--license") + 1], fpm.LICENSE_ARG)


if __name__ == "__main__":
    unittest.main()
