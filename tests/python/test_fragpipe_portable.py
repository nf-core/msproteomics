"""bin/fragpipe_portable.py: FRAGPIPEPLAN's files with every staged input back at its original path.

Run: python3 -m unittest discover -s tests/python
"""

import importlib.util
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("fragpipe_portable", REPO / "bin/fragpipe_portable.py")
fp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fp)

TASK = "/w/ab/cdef"
TOOLS = "/fp/tools"
ROWS = [
    ("lcms", "raw_files/plex1/a.mzML", "s3://bucket/run/a.mzML"),
    ("lcms", "raw_files/plex1/b.mzML", "s3://bucket/run/b.mzML"),
    ("file", "raw_files/plex1/plex1_annotation.txt", "s3://bucket/run/plex1_annotation.txt"),
    ("file", "db.fasta", "/data/db.fasta"),
]


def paths(rows=ROWS):
    return fp.PathMap(TASK, TOOLS, rows)


class PathMapTest(unittest.TestCase):

    def test_inputs_get_their_original_paths(self):
        text = f"database.db-path={TASK}/db.fasta\n{TASK}/raw_files/plex1/a.mzML\tplex1\t\tDDA\n"
        self.assertEqual(paths().apply(text), "database.db-path=/data/db.fasta\ns3://bucket/run/a.mzML\tplex1\t\tDDA\n")

    def test_files_fragpipe_writes_next_to_an_input_go_to_its_original_dir(self):
        self.assertEqual(paths().apply(f"{TASK}/raw_files/plex1/a_uncalibrated.mzML"), "s3://bucket/run/a_uncalibrated.mzML")

    def test_a_file_in_a_staged_dir_with_one_origin_goes_to_that_dir(self):
        self.assertEqual(paths().apply(f"{TASK}/raw_files/plex1/diaumpire_se.log"), "s3://bucket/run/diaumpire_se.log")

    def test_files_in_a_staged_dir_with_two_origins_go_by_stem(self):
        rows = [("lcms", "raw_files/a.mzML", "/x/a.mzML"), ("lcms", "raw_files/b.mzML", "/y/b.mzML")]
        self.assertEqual(paths(rows).apply(f"{TASK}/raw_files/b_uncalibrated.mzML"), "/y/b_uncalibrated.mzML")

    def test_a_name_no_input_stem_starts_in_a_dir_with_two_origins_is_an_error(self):
        rows = [("lcms", "raw_files/a.mzML", "/x/a.mzML"), ("lcms", "raw_files/b.mzML", "/y/b.mzML")]
        with self.assertRaisesRegex(fp.PortableError, "no staged input's stem starts the name"):
            paths(rows).apply(f"{TASK}/raw_files/diaumpire_se.log")

    def test_a_staged_dir_on_its_own_is_its_files_original_dir(self):
        self.assertEqual(paths().apply(f"--dir {TASK}/raw_files/plex1 "), "--dir s3://bucket/run ")

    def test_a_staged_dir_with_files_from_two_dirs_is_an_error(self):
        rows = [("lcms", "raw_files/a.mzML", "/x/a.mzML"), ("lcms", "raw_files/b.mzML", "/y/b.mzML")]
        with self.assertRaisesRegex(fp.PortableError, "named on its own"):
            paths(rows).apply(f"raw_file_location = {TASK}/raw_files\n")

    def test_work_and_tools_dirs_become_placeholders(self):
        text = f"workdir={TASK}/fp_work\nfragpipe-config.tools-folder={TOOLS}\nx={TOOLS}/diann/diann\n"
        self.assertEqual(paths().apply(text), "workdir=@WORKDIR@\nfragpipe-config.tools-folder=@TOOLS_DIR@\nx=@TOOLS_DIR@/diann/diann\n")

    def test_a_task_path_left_is_an_error(self):
        with self.assertRaisesRegex(fp.PortableError, "a task path is left"):
            paths().check("x", f"tmp={TASK}/tmp/x")


class RunLogTest(unittest.TestCase):

    def run_main(self, name, text):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / name
            f.write_text(text)
            m = Path(d) / "map.tsv"
            m.write_text("".join("\t".join(r) + "\n" for r in ROWS))
            out = Path(d) / "out"
            code = fp.main(["--task-dir", TASK, "--tools-dir", TOOLS, "--map", str(m), "--outdir", str(out), str(f)])
            return code, (out / name).read_text() if code == 0 else None

    def test_run_log_gets_original_paths_and_the_licence_redacted(self):
        log = (f"IonQuant [Work dir: {TASK}/fp_work]\n"
               f"java -jar {TOOLS}/IonQuant.jar --license /secret/license.dat --specdir {TASK}/raw_files/plex1 {TASK}/raw_files/plex1/a.mzML\n")
        code, result = self.run_main("log_2026-10-12.txt", log)
        self.assertEqual(code, 0)
        self.assertEqual(result, "IonQuant [Work dir: @WORKDIR@]\n"
                                 "java -jar @TOOLS_DIR@/IonQuant.jar --license {LICENSE} --specdir s3://bucket/run s3://bucket/run/a.mzML\n")

    def test_the_fasta_sample_header_of_the_run_log_gets_the_original_path(self):
        code, result = self.run_main("log_2026-10-12.txt", f"~~~~~~Sample of {TASK}/db.fasta~~~~~~~\n")
        self.assertEqual((code, result), (0, "~~~~~~Sample of /data/db.fasta~~~~~~~\n"))

    def test_licence_material_left_in_another_file_is_an_error(self):
        code, _result = self.run_main("fragpipe.workflow", "x=/secret/license.dat\n")
        self.assertEqual(code, 1)


if __name__ == "__main__":
    unittest.main()
