"""Unit tests for bin/golden_argv_compare.py (run: python3 -m pytest tests/python/test_golden_argv_compare.py)."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

import golden_argv_compare as gac  # noqa: E402

FP = "/fragpipe_bin/fragpipe-24.0/fragpipe-24.0"
JAR = f"{FP}/tools/MSFragger-4.4.1/MSFragger-4.4.1.jar"

RECORDED = ["java", "-jar", "-Dfile.encoding=UTF-8", "-Xmx100G", JAR, "{RUN}/work/split_peptide_index_tempdir/fragger.params",
            "--license", "{LICENSE}", "--split1", "{RUN}/input/run_A01.mzML"]
FP_DIRS = ["{RUN}/work/split_peptide_index_tempdir", "{RUN}/input"]
TASK_DIRS = ["/work/ab/cd"]
MODULE = ["java", "-Duser.home=/work/ab/cd", "-jar", "-Dfile.encoding=UTF-8", "-Xmx64G", JAR, "/work/ab/cd/calibrate_fragger.params",
          "--license", f"{FP}/license.dat", "--split1", "run_B01.mzML"]


class TestNormalise(unittest.TestCase):
    def test_module_command_matches_recorded_command_of_another_run(self):
        recorded = [t.replace("A01", "B01") for t in RECORDED]
        r = gac.compare(recorded, MODULE, FP, "run_B01", "run_B01", FP_DIRS, TASK_DIRS)
        self.assertEqual(r["status"], "MATCH", r)

    def test_missing_license_is_a_diff(self):
        module = MODULE[:7] + MODULE[9:]
        r = gac.compare(RECORDED, module, FP, "run_A01", "run_B01", FP_DIRS, TASK_DIRS)
        self.assertEqual(r["status"], "DIFF")

    def test_extra_native_library_flag_is_a_diff(self):
        module = MODULE[:4] + ["-Dlibs.thermo.dir=/x/ext/thermo"] + MODULE[4:]
        r = gac.compare(RECORDED, module, FP, "run_A01", "run_B01", FP_DIRS, TASK_DIRS)
        self.assertEqual(r["status"], "DIFF")

    def test_wrong_option_order_is_a_diff(self):
        module = MODULE[:2] + MODULE[3:5] + [MODULE[2]] + MODULE[5:]
        r = gac.compare(RECORDED, module, FP, "run_A01", "run_B01", FP_DIRS, TASK_DIRS)
        self.assertEqual(r["status"], "DIFF")

    def test_fragpipe_per_data_type_params_names_are_params(self):
        # FragPipe saves the DDA params as fragger_dda.params when the run has DIA files (CmdMsfragger.java:444)
        for name in ("fragger_dda.params", "fragger_dia.params", "fragger_gpfdia.params", "fragger_dda_plus.params"):
            with self.subTest(name=name):
                recorded = [t.replace("A01", "B01").replace("/fragger.params", "/" + name) for t in RECORDED]
                r = gac.compare(recorded, MODULE, FP, "run_B01", "run_B01", FP_DIRS, TASK_DIRS)
                self.assertEqual(r["status"], "MATCH", r)

    def test_other_params_file_is_not_params(self):
        recorded = [t.replace("A01", "B01").replace("/fragger.params", "/fragger_other.params") for t in RECORDED]
        r = gac.compare(recorded, MODULE, FP, "run_B01", "run_B01", FP_DIRS, TASK_DIRS)
        self.assertEqual(r["status"], "DIFF")

    def test_file_outside_the_task_dirs_is_a_diff(self):
        module = [t.replace("/work/ab/cd/calibrate_fragger.params", "/work/ef/gh/calibrate_fragger.params") for t in MODULE]
        r = gac.compare(RECORDED, module, FP, "run_A01", "run_B01", FP_DIRS, TASK_DIRS)
        self.assertEqual(r["status"], "DIFF")

    def test_spectra_in_another_directory_is_a_diff(self):
        module = MODULE[:-1] + ["/data/other/run_B01.mzML"]
        r = gac.compare(RECORDED, module, FP, "run_A01", "run_B01", FP_DIRS, TASK_DIRS)
        self.assertEqual(r["status"], "DIFF")

    def test_recorded_file_outside_the_fixture_dirs_is_a_diff(self):
        recorded = RECORDED[:-1] + ["{RUN}/elsewhere/run_A01.mzML"]
        r = gac.compare(recorded, MODULE, FP, "run_A01", "run_B01", FP_DIRS, TASK_DIRS)
        self.assertEqual(r["status"], "DIFF")



if __name__ == "__main__":
    unittest.main()
