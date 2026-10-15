"""Tests for bin/fragpipe_second_pass.py (FragPipe 24.0 finalizer, FragpipeRun.java:589-700)."""

import importlib.util
import os
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "fragpipe_second_pass.py"
_spec = importlib.util.spec_from_file_location("fragpipe_second_pass", SCRIPT)
sp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sp)

GOLDEN = os.environ.get("GOLDEN_ROOT")


class TestJavaFloat:
    # Float.toString(Float.parseFloat(s)) as JDK 17 prints them
    @pytest.mark.parametrize("text,expected", [
        ("10.000000", "10.0"), ("0.000000", "0.0"), ("7", "7.0"), ("9.50", "9.5"), ("0.01", "0.01"),
        ("12.5", "12.5"), ("20", "20.0"), ("0.0001", "1.0E-4"), ("12345678", "1.2345678E7"), ("0.001", "0.001"),
    ])
    def test_java_float_to_string(self, text, expected):
        assert sp.java_float_to_string(text) == expected


class TestConsole:
    def test_first_match_and_zero_width_space_lines_dropped(self):
        text = "x New fragment_mass_tolerance = 1 PPM​\nNew fragment_mass_tolerance = 7 PPM\nNew fragment_mass_tolerance = 9.5 PPM\n"
        assert sp.console_values(text)["fragment_mass_tolerance"] == "7.0"

    def test_tolerance_needs_ppm(self):
        assert sp.console_values("New fragment_mass_tolerance = 0.02 Da\n")["fragment_mass_tolerance"] is None


class TestWorkflow:
    def test_rules(self):
        wf = "\n".join([
            "# FragPipe workflow", "  database.db-path=/data/db.fasta  ", "database.decoy-tag=rev_",
            "msfragger.calibrate_mass=2", "msfragger.fragment_mass_tolerance=20", "msfragger.fragment_mass_units=1",
            "msfragger.use_topN_peaks=150", "msfragger.minimum_ratio=0.01", "msfragger.intensity_transform=0",
            "msfragger.remove_precursor_peak=1", "msfragger.max_fragment_charge=2", "msfragger.precursor_true_tolerance=20",
            "msfragger.precursor_true_units=1", "peptide-prophet.run-peptide-prophet=true", "percolator.run-percolator=false",
            "phi-report.print-decoys=true", "tab-run.sub_mzml_prob_threshold=0.5", "tab-run.write_sub_mzml=true",
            "workdir=/out/fragpipe", "",
        ])
        console = "\n".join(["progress line​", "New fragment_mass_tolerance = 12 PPM", "New use_topN_peaks = 100",
                             "New minimum_ratio = 0.00", "New intensity_transform = 1", "New remove_precursor_peak = 0",
                             "New fragment_mass_tolerance = 30 PPM", ""])
        assert sp.second_pass_workflow(wf, console, "/fp/fragpipe-24.0/databases") == "\n".join([
            "# FragPipe workflow", "database.db-path=/fp/fragpipe-24.0/databases", "database.decoy-tag=",
            "msfragger.calibrate_mass=0", "msfragger.fragment_mass_tolerance=12.0", "msfragger.fragment_mass_units=1",
            "msfragger.use_topN_peaks=100", "msfragger.minimum_ratio=0.0", "msfragger.intensity_transform=1",
            "msfragger.remove_precursor_peak=0", "msfragger.max_fragment_charge=2", "msfragger.precursor_true_tolerance=20",
            "msfragger.precursor_true_units=1", "peptide-prophet.run-peptide-prophet=false", "percolator.run-percolator=true",
            "phi-report.print-decoys=false", "tab-run.sub_mzml_prob_threshold=0.5", "tab-run.write_sub_mzml=false",
            "workdir=", "",
        ])

    def test_units_stay_without_a_tolerance(self):
        wf = "msfragger.fragment_mass_tolerance=0.02\nmsfragger.fragment_mass_units=0\n"
        assert sp.second_pass_workflow(wf, "New fragment_mass_tolerance = 0.02 Da\n", "/d") == wf

    def test_crlf_lines(self):
        assert sp.second_pass_workflow("a=1\r\nb=2", "", "/d") == "a=1\nb=2\n"


class TestManifest:
    def test_sub_mzml_paths_and_java_split(self):
        text = "/in/A01.mzML\tA\t1\tDDA\n/in/run.b.raw\tB\t\tDDA"
        assert sp.second_pass_manifest(text, "/out/./fragpipe") == \
            "/out/fragpipe/A01_sub.mzML\tA\t1\tDDA\n/out/fragpipe/run.b_sub.mzML\tB\t\tDDA\n"

    # wd.resolve(name).toAbsolutePath().normalize() (FragpipeRun.java:694) on the path; a URI keeps its scheme and bucket
    @pytest.mark.parametrize("workdir,expected", [
        ("s3://b/o/fragpipe", "s3://b/o/fragpipe/a_sub.mzML"),
        ("s3://b/o/./x/../fragpipe/", "s3://b/o/fragpipe/a_sub.mzML"),
        ("/data/run/../out/fragpipe", "/data/out/fragpipe/a_sub.mzML"),
    ])
    def test_uri_and_local_work_dirs(self, workdir, expected):
        assert sp.second_pass_manifest("/in/a.mzML\tA\t1\tDDA", workdir) == expected + "\tA\t1\tDDA\n"

    def test_relative_work_dir_fails(self):
        with pytest.raises(ValueError, match="absolute path or a URI"):
            sp.second_pass_manifest("/in/a.mzML\tA\t1\tDDA", "out/fragpipe")

    def test_trailing_empty_fields_dropped(self):
        assert sp.second_pass_manifest("/in/a.mzML\t\t\tDDA\t\t\n", "/w") == "/w/a_sub.mzML\t\t\tDDA\n"


@pytest.mark.skipif(not GOLDEN, reason="GOLDEN_ROOT not set")
class TestGolden:
    """Byte-equal to FragPipe 24.0's files of golden set PXD001468/Basic-Search_submzml."""

    def test_byte_equal(self, tmp_path):
        import json
        g = Path(GOLDEN) / "PXD001468" / "Basic-Search_submzml"
        roots = json.loads((g / "manifest.json").read_text())["roots"]
        sp.main(["--workflow", str(g / "files/work/fragpipe.workflow"), "--manifest", str(g / "files/work/fragpipe-files.fp-manifest"),
                 "--log", str(g / "fragpipe_log.redacted.txt"), "--workdir", roots["{WORKDIR}"],
                 "--databases-dir", "/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/databases", "--outdir", str(tmp_path)])
        for name in ("fragpipe-second-pass.workflow", "fragpipe-files-second-pass.fp-manifest"):
            assert (tmp_path / name).read_bytes() == (g / "files/work" / name).read_bytes(), name
