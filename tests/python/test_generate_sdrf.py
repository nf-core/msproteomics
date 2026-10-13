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
