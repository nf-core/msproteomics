"""Tests for bin/generate_sdrf.py: the bookkeeping SDRF uses FragPipe's organism guess."""

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN))
_spec = importlib.util.spec_from_file_location("generate_sdrf", BIN / "generate_sdrf.py")
gs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gs)
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "iproteomes"


class TestOrganism(unittest.TestCase):
    def organism_column(self, organism: str, database: Path) -> list[str]:
        with tempfile.TemporaryDirectory() as tmp:
            sheet = Path(tmp) / "s.csv"
            sheet.write_text("sample,spectra,condition\nA,a.mzML,c1\nB,b.mzML,c2\n")
            out = Path(tmp) / "out.sdrf.tsv"
            gs.generate_sdrf(str(sheet), {"organism": organism}, str(out), str(database))
            with open(out, newline="") as fh:
                return [r["characteristics[organism]"] for r in csv.DictReader(fh, delimiter="\t")]

    def test_given_organism_is_written_as_is(self):
        self.assertEqual(self.organism_column("Mus musculus", FIXTURES / "organisms.fasta"), ["Mus musculus"] * 2)

    def test_empty_organism_takes_the_most_frequent_os_of_the_database(self):
        self.assertEqual(self.organism_column("", FIXTURES / "organisms.fasta"), ["Saccharomyces cerevisiae"] * 2)

    def test_database_without_os_gives_not_available(self):
        self.assertEqual(self.organism_column("", FIXTURES / "test.fasta"), ["not available"] * 2)

    def test_guess_is_the_fragpipe_sdrf_implementation(self):
        self.assertIs(gs.fragpipe_sdrf.guess_organism, sys.modules["fragpipe_sdrf"].guess_organism)


if __name__ == "__main__":
    unittest.main()


class TestDataFile(unittest.TestCase):
    """comment[data file] is the decompressed run file (what DIA-NN / FragPipe and qpx name the run)."""

    # assets/schema_input.json spectra pattern: format and optional archive suffix
    FORMATS = ["raw", "RAW", "mzML", "mzXML", "mgf", "mzBIN", "d", "dia"]
    ARCHIVES = ["", ".gz", ".tar", ".tar.gz", ".zip"]

    def columns(self, spectra: str) -> dict:
        with tempfile.TemporaryDirectory() as tmp:
            sheet = Path(tmp) / "s.csv"
            sheet.write_text(f"sample,spectra,condition\nA,{spectra},c1\n")
            out = Path(tmp) / "out.sdrf.tsv"
            gs.generate_sdrf(str(sheet), {"organism": "Homo sapiens"}, str(out), str(FIXTURES / "test.fasta"))
            with open(out, newline="") as fh:
                return next(csv.DictReader(fh, delimiter="\t"))

    def test_every_schema_combination_names_the_decompressed_file(self):
        import qpx_export_lib

        for fmt in self.FORMATS:
            for archive in self.ARCHIVES:
                with self.subTest(fmt=fmt, archive=archive):
                    row = self.columns(f"s3://bucket/run/sample.{fmt}{archive}")
                    self.assertEqual(row["comment[data file]"], f"sample.{fmt}")
                    self.assertEqual(row["comment[file uri]"], f"s3://bucket/run/sample.{fmt}{archive}")
                    # the qpx export's run of the SDRF row is the run the search engine reports
                    expected = "sample" if fmt.lower() in ("raw", "mzml", "mzxml", "mgf", "d", "dia") else f"sample.{fmt}"
                    self.assertEqual(qpx_export_lib.run_file_stem(row["comment[data file]"]), expected)
