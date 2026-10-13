"""Unit tests for bin/golden_compare.py (run: python3 -m unittest discover tests/python)."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

import golden_compare as gc  # noqa: E402

UUID = "2c5b7a2f-bf30-4a8c-a1f4-d1a77bcb6527"
UUID2 = "397f4c49-2094-4394-ad02-2db4b8f98d2c"


class Base(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.t = Path(self._tmp.name)
        self.exp = gc.Side(["/runs/a"], ["/runs/a/jtmp"])
        self.act = gc.Side(["/runs/b"], ["/runs/b/jtmp"])

    def tearDown(self):
        self._tmp.cleanup()

    def pair(self, name, a: bytes, b: bytes, profile: str = "dda", ionquant_rewritten=None) -> gc.Result:
        (self.t / "e").mkdir(exist_ok=True)
        (self.t / "a").mkdir(exist_ok=True)
        pe, pa = self.t / "e" / Path(name).name, self.t / "a" / Path(name).name
        pe.write_bytes(a)
        pa.write_bytes(b)
        return gc.compare_files(pe, pa, name, self.exp, self.act, profile, ionquant_rewritten)


class TestRuleSelection(unittest.TestCase):
    def test_rules_by_file_name(self):
        cases = {
            "log_2026-10-04_23-41-38.txt": gc.SKIP,
            "x.fasta.1.pepindex": gc.PEPINDEX,
            "combined.prot.xml": gc.PROTXML,
            "a.pepXML": gc.PEPXML,
            "interact-a.pep.xml": gc.PEPXML,
            "filter.log": gc.PHILO_FILTER_LOG,
            "fragger.params": gc.PATH_ONLY,
            "msbooster_params.txt": gc.PATH_ONLY,
            "fragpipe.workflow": gc.WORKFLOW_FILE,
            "fragpipe.job": gc.PATH_ONLY,
            "fragpipe-files.fp-manifest": gc.PATH_ONLY,
            "filelist_ionquant.txt": gc.PATH_ONLY,
            "experiment_annotation.tsv": gc.PATH_ONLY,
            "psm.tsv": gc.PHILO_PSM,
            "peptide.tsv": gc.PHILO_COUNTS,
            "ion.tsv": gc.PHILO_COUNTS,
            "protein.tsv": gc.PHILO_COUNTS,
            "library.tsv": gc.RAW,
            "a.pin": gc.RAW,
            "msstats.csv": gc.PHILO_RUN,
            "protein.fas": gc.RAW,
            "a_uncalibrated.mzBIN": gc.RAW,
        }
        for name, rule in cases.items():
            with self.subTest(name=name):
                self.assertEqual(gc.rule_for(name), rule)


class TestRaw(Base):
    def test_identical_bytes_match(self):
        self.assertEqual(self.pair("library.tsv", b"a\tb\n1\t2\n", b"a\tb\n1\t2\n").status, "MATCH")

    def test_difference_reports_offset_and_line(self):
        r = self.pair("library.tsv", b"a\tb\n1\t2\n", b"a\tb\n1\t3\n")
        self.assertEqual(r.status, "DIFF")
        self.assertIn("byte 6, line 2", r.detail)

    def test_paths_in_tsv_are_not_normalized(self):
        self.assertEqual(self.pair("library.tsv", b"/runs/a/x\n", b"/runs/b/x\n").status, "DIFF")


class TestPepXML(Base):
    def test_path_and_iso_timestamps_normalized(self):
        a = b'<msms_pipeline_analysis date="2026-10-04T23:38:10" summary_xml="/runs/a/input/x.pepXML">\n<analysis_summary analysis="Percolator" time="2026-10-04T23:41:18">\n'
        b = b'<msms_pipeline_analysis date="2026-10-05T01:00:00" summary_xml="/runs/b/input/x.pepXML">\n<analysis_summary analysis="Percolator" time="2026-10-05T01:02:03">\n'
        self.assertEqual(self.pair("interact-x.pep.xml", a, b).status, "MATCH")

    def test_value_difference_still_detected(self):
        a = b'<parameter name="intensity_transform" value="1"/>\n'
        b = b'<parameter name="intensity_transform" value="0"/>\n'
        r = self.pair("x.pepXML", a, b)
        self.assertEqual(r.status, "DIFF")
        self.assertIn("line 1", r.detail)

    def test_non_iso_time_is_not_normalized(self):
        self.assertEqual(self.pair("x.pepXML", b'time="12:00:00"\n', b'time="12:00:01"\n').status, "DIFF")

    def test_prefix_only_matches_at_path_boundary(self):
        self.assertEqual(self.pair("x.pepXML", b'f="/runs/a/x"\n', b'f="/runs/bb/x"\n').status, "DIFF")


class TestPepXMLTmp(Base):
    def test_iprophet_summary_xml_tmp_uuid_dir_normalized(self):
        a = f'<msms_pipeline_analysis date="2026-10-06T20:22:51" summary_xml="/runs/a/jtmp/{UUID}/combined.pep.xml">\n'
        b = f'<msms_pipeline_analysis date="2026-10-07T00:10:13" summary_xml="/runs/b/jtmp/{UUID2}/combined.pep.xml">\n'
        self.assertEqual(self.pair("combined.pep.xml", a.encode(), b.encode()).status, "MATCH")

    def test_tmp_dir_outside_the_tmp_roots_is_a_diff(self):
        a = f'<x summary_xml="/runs/a/jtmp/{UUID}/combined.pep.xml">\n'
        b = f'<x summary_xml="/elsewhere/{UUID2}/combined.pep.xml">\n'
        self.assertEqual(self.pair("combined.pep.xml", a.encode(), b.encode()).status, "DIFF")


class TestFpMetaConfig(Base):
    def test_db_optimizer_config_run_paths_are_normalized(self):
        a = b"projectDir = /runs/a/work\nfastaFile = /runs/a/db/db.fasta\noutdir = /runs/a/work\nqvalue = 0.01\n"
        b = b"projectDir = /runs/b/work\nfastaFile = /runs/b/db/db.fasta\noutdir = /runs/b/work\nqvalue = 0.01\n"
        self.exp = gc.Side(["/runs/a/work", "/runs/a/db"])
        self.act = gc.Side(["/runs/b/work", "/runs/b/db"])
        self.assertEqual(self.pair("DbOptimizerConfig.yaml", a, b).status, "MATCH")

    def test_db_optimizer_config_value_difference_is_a_diff(self):
        self.exp = gc.Side(["/runs/a"])
        self.act = gc.Side(["/runs/b"])
        self.assertEqual(self.pair("DbOptimizerConfig.yaml", b"qvalue = 0.01\n", b"qvalue = 0.02\n").status, "DIFF")


class TestFpopTables(Base):
    HEAD = "Group Name\tProtein\tPeptide Sequence\tStart\tFPOP Mods\tOther Mods\t\n"

    def test_peptide_list_order_is_not_a_difference(self):
        a = self.HEAD + "G1\tP1\tAAK, BBK, CCK\t3\t7W(15.9949)\t\t\n"
        b = self.HEAD + "G1\tP1\tCCK, AAK, BBK\t3\t7W(15.9949)\t\t\n"
        self.assertEqual(self.pair("FPOP_peptides.tsv", a.encode(), b.encode()).status, "MATCH")

    def test_different_peptides_in_the_list_are_a_diff(self):
        a = self.HEAD + "G1\tP1\tAAK, BBK\t3\t7W(15.9949)\t\t\n"
        b = self.HEAD + "G1\tP1\tAAK, CCK\t3\t7W(15.9949)\t\t\n"
        self.assertEqual(self.pair("FPOP_sites.tsv", a.encode(), b.encode()).status, "DIFF")

    def test_order_of_the_other_columns_is_not_normalized(self):
        a = self.HEAD + "G1\tP1\tAAK\t3\t7W(15.9949), 6L(15.9949)\t\t\n"
        b = self.HEAD + "G1\tP1\tAAK\t3\t6L(15.9949), 7W(15.9949)\t\t\n"
        self.assertEqual(self.pair("FPOP_peptides.tsv", a.encode(), b.encode()).status, "DIFF")

    def test_row_count_and_header_are_compared(self):
        a = self.HEAD + "G1\tP1\tAAK\t3\t\t\t\n"
        self.assertEqual(self.pair("FPOP_peptides.tsv", a.encode(), (a + a.split("\n", 1)[1]).encode()).status, "DIFF")
        self.assertEqual(self.pair("FPOP_peptides.tsv", a.encode(), a.replace("Protein", "Prot").encode()).status, "DIFF")


class TestProtXML(Base):
    def test_tmp_uuid_path_and_time_normalized(self):
        a = f'<?xml-stylesheet href="/runs/a/jtmp/{UUID}/combined.prot.xsl"?>\n<program_details time="2026-10-04T23:41:20"/>\n<h source_files="/runs/a/work/C02/i.pep.xml"/>\n'
        b = f'<?xml-stylesheet href="/runs/b/jtmp/{UUID2}/combined.prot.xsl"?>\n<program_details time="2026-10-04T23:41:14"/>\n<h source_files="/runs/b/work/C02/i.pep.xml"/>\n'
        self.assertEqual(self.pair("combined.prot.xml", a.encode(), b.encode()).status, "MATCH")

    def test_subdirectory_difference_is_a_diff(self):
        a = b'<h source_files="/runs/a/C02/i.pep.xml"/>\n'
        b = b'<h source_files="/runs/b/i.pep.xml"/>\n'
        self.assertEqual(self.pair("combined.prot.xml", a, b).status, "DIFF")


class TestPhilosopherLogAndPepindex(Base):
    def test_filter_log_time_stripped(self):
        a = b'time="23:41:21" level=info msg="Parsing /runs/a/work/C02/i.pep.xml"\n'
        b = b'time="23:41:30" level=info msg="Parsing /runs/b/work/C02/i.pep.xml"\n'
        self.assertEqual(self.pair("filter.log", a, b).status, "MATCH")

    def test_filter_log_count_difference_detected(self):
        a = b'time="23:41:21" level=info msg="x" decoy=118 target=1578\n'
        b = b'time="23:41:21" level=info msg="x" decoy=119 target=1578\n'
        self.assertEqual(self.pair("filter.log", a, b).status, "DIFF")

    def test_pepindex_header_and_path_ignored_body_compared(self):
        a = b"HDR1\x00\x00\x00\x20/runs/a/db/x.fasta" + b"BODY"
        b = b"HDR2\x00\x00\x00\x21/runs/b/db/x.fasta" + b"BODY"
        self.assertEqual(self.pair("x.fasta.1.pepindex", a, b).status, "MATCH")
        self.assertEqual(self.pair("x.fasta.2.pepindex", a, b[:-1] + b"Z").status, "DIFF")


class TestDirs(Base):
    def test_missing_extra_and_skipped_files(self):
        e, a = self.t / "E", self.t / "A"
        for d in (e, a):
            (d / "C02").mkdir(parents=True)
            (d / "C02" / "library.tsv").write_text("x\n")
        (e / "only_expected.tsv").write_text("1\n")
        (a / "only_actual.tsv").write_text("1\n")
        (e / "log_2026-10-04_23-41-38.txt").write_text("t\n")
        res = {r.path: r.status for r in gc.compare_dirs(e, a, self.exp, self.act)}
        self.assertEqual(res, {"C02/library.tsv": "MATCH", "only_expected.tsv": "MISSING", "only_actual.tsv": "EXTRA", "log_2026-10-04_23-41-38.txt": "SKIP"})
        self.assertTrue(gc.failed(gc.compare_dirs(e, a, self.exp, self.act)))

    def test_cli_exit_codes(self):
        e, a = self.t / "E", self.t / "A"
        e.mkdir()
        a.mkdir()
        (e / "x.tsv").write_text("1\n")
        (a / "x.tsv").write_text("1\n")
        self.assertEqual(gc.main(["dirs", str(e), str(a)]), 0)
        (a / "x.tsv").write_text("2\n")
        self.assertEqual(gc.main(["dirs", str(e), str(a), "--json", str(self.t / "r.json")]), 1)
        self.assertTrue((self.t / "r.json").is_file())

    def test_pairs_mode_missing_actual_is_failure(self):
        (self.t / "x.tsv").write_text("1\n")
        self.assertEqual(gc.main(["pairs", "--pair", "x.tsv", str(self.t / "x.tsv"), str(self.t / "nope.tsv")]), 1)

    def test_two_empty_trees_fail(self):
        e, a = self.t / "E", self.t / "A"
        e.mkdir()
        a.mkdir()
        res = gc.compare_dirs(e, a, self.exp, self.act)
        self.assertEqual([r.status for r in res], ["EMPTY"])
        self.assertEqual(gc.main(["dirs", str(e), str(a)]), 1)

    def test_only_skipped_files_fail(self):
        e, a = self.t / "E", self.t / "A"
        for d in (e, a):
            d.mkdir()
            (d / "log_2026-10-04_23-41-38.txt").write_text("x\n")
        self.assertTrue(gc.failed(gc.compare_dirs(e, a, self.exp, self.act)))

    def test_question_mark_dir_skipped(self):
        e, a = self.t / "E", self.t / "A"
        for d in (e, a):
            (d / "?" / ".java").mkdir(parents=True)
            (d / "x.tsv").write_text("1\n")
        (e / "?" / ".java" / "font.cache").write_text("a")
        res = {r.path: r.status for r in gc.compare_dirs(e, a, self.exp, self.act)}
        self.assertEqual(res["?/.java/font.cache"], "SKIP")
        self.assertFalse(gc.failed(gc.compare_dirs(e, a, self.exp, self.act)))


class TestDiatracerMzml(Base):
    HEAD = b'<sourceFile location="%s">\n<run startTimeStamp="%s">\n<binary>AAAA</binary>\n'
    TAIL = b'<offset idRef="frame=1">%d</offset>\n</index>\n<indexListOffset>%d</indexListOffset>\n'

    def make(self, run, stamp, off, body=b"<binary>AAAA</binary>\n"):
        return (self.HEAD % (run, stamp)).replace(b"<binary>AAAA</binary>\n", body) + self.TAIL % (off, off + 7)

    def test_run_path_time_and_offsets_are_accepted(self):
        a = self.make(b"/runs/a", b"2026-10-07T05:13:29.273Z", 1000)
        b = self.make(b"/runs/b", b"2026-10-07T18:20:08.118Z", 1125)
        self.assertEqual(self.pair("input/x_diatracer.mzML", a, b).status, "MATCH")

    def test_spectrum_data_is_exact(self):
        a = self.make(b"/runs/a/input", b"2026-10-07T05:13:29.273Z", 1000)
        b = self.make(b"/runs/b/input", b"2026-10-07T05:13:29.273Z", 1000, body=b"<binary>AAAB</binary>\n")
        self.assertEqual(self.pair("input/x_diatracer.mzML", a, b).status, "DIFF")

    def test_line_count_is_exact(self):
        a = self.make(b"/runs/a/input", b"t", 1000)
        self.assertEqual(self.pair("input/x_diatracer.mzML", a, a + b"<extra/>\n").status, "DIFF")

    def test_write_sub_mzml_output_has_the_same_rule(self):
        a = self.make(b"/runs/a", b"t", 1000)
        b = self.make(b"/runs/b", b"t", 1010)
        self.assertEqual(self.pair("work/x_sub.mzML", a, b).status, "MATCH")
        c = self.make(b"/runs/b", b"t", 1010, body=b"<binary>AAAB</binary>\n")
        self.assertEqual(self.pair("work/x_sub.mzML", a, c).status, "DIFF")


class TestDiaProfile(Base):
    """Rules from the determinism doc DIA sections; each one only in the dia profile."""

    def test_dia_rules_by_path(self):
        cases = {
            "MSBooster/spectraRT.predicted.bin": gc.SKIP,
            "easypqp_files/easypqp_rt_alignment_x.pdf": gc.PDF,
            "dia-quant-output/report.log.txt": gc.DIANN_LOG,
            "filter.log": gc.PHILO_FILTER_LOG,
            "MSBooster/spectraRT.tsv": gc.SORT_ROWS,
            "MSBooster/spectraRT_full.tsv": gc.SORT_ROWS,
            "psm.tsv": gc.PHILO_PSM,
            "peptide.tsv": gc.PHILO_FULL,
            "ion.tsv": gc.PHILO_FULL,
            "protein.tsv": gc.PHILO_COUNTS,
            "msbooster_params.txt": gc.PATH_ONLY,
            "easypqp_files/cmds.txt": gc.EASYPQP_CMDS,
            "library.tsv": gc.RAW,
            "dia-quant-output/report.tsv": gc.PATH_ONLY,
            "dia-quant-output/report.pg_matrix.tsv": gc.PATH_ONLY,
            "dia-quant-output/report.pr_matrix.tsv": gc.PATH_ONLY,
            "dia-quant-output/report.gg_matrix.tsv": gc.PATH_ONLY,
            "dia-quant-output/report.unique_genes_matrix.tsv": gc.PATH_ONLY,
            "dia-quant-output/report.stats.tsv": gc.PATH_ONLY,
            "dia-quant-output/report-tsv.speclib": gc.RAW,
            "dia-quant-output/msstats.csv": gc.RAW,
            "report.tsv": gc.RAW,
            "dia-quant-output/abundance_peptide_MS2quant_Norm.tsv": gc.FRAGREPORTER_TABLE,
            "tmt-report/abundance_protein_MD.tsv": gc.RAW,
        }
        for rel, rule in cases.items():
            with self.subTest(rel=rel):
                self.assertEqual(gc.rule_for(rel, "dia"), rule)

    def test_dda_profile_keeps_dia_files_strict(self):
        for rel in ("MSBooster/spectraRT.tsv", "x.pdf", "dia-quant-output/report.log.txt", "MSBooster/spectraRT.predicted.bin",
                    "dia-quant-output/report.tsv", "dia-quant-output/abundance_peptide_MS1quant_None.tsv"):
            with self.subTest(rel=rel):
                self.assertEqual(gc.rule_for(rel, "dda"), gc.RAW)
        # IonQuant rewrites peptide/ion in DDA: only the spectral-count rule (doc "Workflow scope")
        self.assertEqual(gc.rule_for("ion.tsv", "dda"), gc.PHILO_COUNTS)

    def test_manifest_decides_peptide_ion_envelope_not_profile(self):
        for profile in ("dda", "dia"):
            self.assertEqual(gc.rule_for("ion.tsv", profile, ionquant_rewritten=False), gc.PHILO_FULL)
            self.assertEqual(gc.rule_for("peptide.tsv", profile, ionquant_rewritten=True), gc.PHILO_COUNTS)

    def test_tmt_pre_ionquant_ion_tsv_within_bounds_passes_but_strict_when_rewritten(self):
        a = (b"Peptide Sequence\tModified Sequence\tPrev AA\tNext AA\tCharge\tM/Z\tExpectation\tProbability\tSpectral Count\tQvalue\n"
             b"AAAK\tAAAK\tR\tS\t2\t500.1\t0.001\t0.99\t3\t0\n")
        b = a.replace(b"0.001", b"0.002")
        self.assertEqual(self.pair("ion.tsv", a, b, "dda", ionquant_rewritten=False).status, "MATCH")
        self.assertEqual(self.pair("ion.tsv", a, b, "dda", ionquant_rewritten=True).status, "DIFF")

    def test_unknown_profile_raises(self):
        with self.assertRaises(ValueError):
            gc.rule_for("psm.tsv", "dial")

    def test_pdf_creation_date(self):
        a = b"%PDF\n/CreationDate (D:20261005191642Z)\n/Producer (matplotlib)\n"
        b = b"%PDF\n/CreationDate (D:20261005191743Z)\n/Producer (matplotlib)\n"
        self.assertEqual(self.pair("easypqp_files/a.pdf", a, b, "dia").status, "MATCH")
        self.assertEqual(self.pair("easypqp_files/a.pdf", a, b.replace(b"matplotlib", b"other"), "dia").status, "DIFF")

    def test_diann_log(self):
        a = b"DIA-NN 1.8\nCurrent date and time: Mon Oct  5 19:16:43 2026\n/runs/a/work/x --cfg\n[0:01] Loading\n[1:20] Done\n"
        b = b"DIA-NN 1.8\nCurrent date and time: Mon Oct  5 19:17:44 2026\n/runs/b/work/x --cfg\n[0:02] Loading\n[1:31] Done\n"
        self.assertEqual(self.pair("dia-quant-output/report.log.txt", a, b, "dia").status, "MATCH")
        self.assertEqual(self.pair("dia-quant-output/report.log.txt", a, b.replace(b"Done", b"Fail"), "dia").status, "DIFF")

    def test_workflow_diann_version_case(self):
        a = b"# DIA-NN version 1.8.2 Beta 8\nworkdir=/runs/a/work\n"
        b = b"# DIA-NN version 1.8.2 beta 8\nworkdir=/runs/b/work\n"
        self.assertEqual(self.pair("fragpipe.workflow", a, b, "dda").status, "MATCH")
        self.assertEqual(self.pair("fragpipe.workflow", a, b.replace(b"1.8.2", b"1.9"), "dda").status, "DIFF")

    def test_diann_log_host_core_count(self):
        a = b"DIA-NN 1.8\nLogical CPU cores: 128\n--threads 64\n"
        b = b"DIA-NN 1.8\nLogical CPU cores: 96\n--threads 64\n"
        self.assertEqual(self.pair("dia-quant-output/report.log.txt", a, b, "dia").status, "MATCH")
        self.assertEqual(self.pair("dia-quant-output/report.log.txt", a, b.replace(b"--threads 64", b"--threads 32"), "dia").status, "DIFF")

    def test_spectrart_row_order_header_fixed(self):
        a = b"peptide\tcharge\nAAA\t2\nBBB\t3\n"
        b = b"peptide\tcharge\nBBB\t3\nAAA\t2\n"
        self.assertEqual(self.pair("MSBooster/spectraRT.tsv", a, b, "dia").status, "MATCH")
        self.assertEqual(self.pair("MSBooster/spectraRT.tsv", a, b"charge\tpeptide\nAAA\t2\nBBB\t3\n", "dia").status, "DIFF")


    def test_diann_report_run_paths(self):
        a = b"File.Name\tRun\tPrecursor.Quantity\n/runs/a/input/r1.d\tr1\t10.5\n"
        b = b"File.Name\tRun\tPrecursor.Quantity\n/runs/b/r1.d\tr1\t10.5\n"
        self.exp = gc.Side(["/runs/a/input", "/runs/a/work"])
        self.assertEqual(self.pair("dia-quant-output/report.tsv", a, b, "dia").status, "MATCH")
        self.assertEqual(self.pair("dia-quant-output/report.stats.tsv", a, b.replace(b"10.5", b"10.6"), "dia").status, "DIFF")

    def test_fragreporter_run_column_order_is_not_a_difference(self):
        # FragReporter orders the run columns by HashMap<absolute path> (DiannPrReport.runInfo): the golden
        # had r2 first, the Tower run r1 first, with the same values per run
        a = b"Peptide\tGene\t/runs/a/input/r2.d\t/runs/a/input/r1.d\nAAK\tG1\t2.0\t1.0\nCCK\tG2\t\t3.0\n"
        b = b"Peptide\tGene\t/runs/b/r1.d\t/runs/b/r2.d\nAAK\tG1\t1.0\t2.0\nCCK\tG2\t3.0\t\n"
        self.exp = gc.Side(["/runs/a/input"])
        rel = "dia-quant-output/abundance_peptide_MS2quant_Norm.tsv"
        self.assertEqual(self.pair(rel, a, b, "dia").status, "MATCH")

    def test_fragreporter_value_swapped_between_runs_is_a_diff(self):
        a = b"Peptide\t/runs/a/r2.d\t/runs/a/r1.d\nAAK\t2.0\t1.0\n"
        b = b"Peptide\t/runs/b/r1.d\t/runs/b/r2.d\nAAK\t2.0\t1.0\n"
        self.assertEqual(self.pair("dia-quant-output/abundance_peptide_MS1quant_None.tsv", a, b, "dia").status, "DIFF")

    def test_fragreporter_without_run_columns_raises(self):
        a = b"Peptide\t/elsewhere/r1.d\nAAK\t1.0\n"
        with self.assertRaises(ValueError):
            self.pair("dia-quant-output/abundance_peptide_MS1quant_None.tsv", a, a, "dia")

    def test_fragreporter_run_columns_not_at_the_end_raise(self):
        a = b"Peptide\t/runs/a/r1.d\tGene\nAAK\t1.0\tG\n"
        b = b"Peptide\t/runs/b/r1.d\tGene\nAAK\t1.0\tG\n"
        with self.assertRaises(ValueError):
            self.pair("dia-quant-output/abundance_peptide_MS1quant_None.tsv", a, b, "dia")

    def test_easypqp_cmds_producer_dirs_are_dropped(self):
        # gen_con_spec_lib.py resolves the staged symlinks to their producer task dirs
        a = (b"python3 easypqp convert --pepxml '['\"'\"'/runs/a/work/interact-r1.pep.xml'\"'\"']' --spectra /runs/a/input/r1.mzML --psms r1.psmpkl\n"
             b"python3 easypqp library --psmtsv /runs/a/work/psm.tsv --out lib.tsv /runs/a/work/filelist_easypqp_library.txt\n")
        b = (b"python3 easypqp convert --pepxml '['\"'\"'/runs/b/t1/r1/interact-r1.pep.xml'\"'\"']' --spectra /runs/b/t2/r1/r1.mzML --psms r1.psmpkl\n"
             b"python3 easypqp library --psmtsv /runs/b/t3/workdir-root/psm.tsv --out lib.tsv /runs/b/t4/filelist_easypqp_library.txt\n")
        self.exp = gc.Side(["/runs/a/work", "/runs/a/input"])
        self.act = gc.Side(["/runs/b/t1", "/runs/b/t2", "/runs/b/t3", "/runs/b/t4"])
        self.assertEqual(self.pair("easypqp_files/cmds.txt", a, b, "dia").status, "MATCH")
        self.assertEqual(self.pair("easypqp_files/cmds.txt", a, b.replace(b"r1.mzML", b"r2.mzML"), "dia").status, "DIFF")
        self.assertEqual(self.pair("easypqp_files/cmds.txt", a, b.replace(b"lib.tsv", b"x.tsv"), "dia").status, "DIFF")


if __name__ == "__main__":
    unittest.main()
