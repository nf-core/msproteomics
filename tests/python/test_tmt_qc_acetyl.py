"""N-terminal acetyl detection of the TMT label check QC: only N-terminal 42.0106
counts as biological blocking; K-acetyl does not."""

import importlib.util
import sys
import types
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "modules/local/tmtlabelcheckanalyze/resources/usr/bin"


def load(name):
    for heavy in ("matplotlib", "matplotlib.pyplot", "numpy", "pandas"):
        try:
            __import__(heavy)
        except ImportError:
            sys.modules.setdefault(heavy, types.ModuleType(heavy))
    if not hasattr(sys.modules["matplotlib"], "use"):
        sys.modules["matplotlib"].use = lambda *a, **k: None
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestPsmAcetyl(unittest.TestCase):
    qc = load("tmt_qc")

    def test_nterm_acetyl_assigned_modifications(self):
        # golden C02/psm.tsv row
        self.assertTrue(self.qc.has_nterm_acetylation("N-term(42.0106),12C(57.0215),8M(15.9949)"))

    def test_lysine_acetyl_is_not_nterm(self):
        self.assertFalse(self.qc.has_nterm_acetylation("N-term(304.2071),5K(42.0106)"))

    def test_philosopher_modified_peptide(self):
        self.assertTrue(self.qc.has_nterm_acetylation("n[43]SNGYEDHM[147]AEDCR"))
        self.assertFalse(self.qc.has_nterm_acetylation("n[305]PEPK[170]R"))


class TestIonQuantAcetyl(unittest.TestCase):
    qc = load("tmt_qc_ionquant")

    def test_nterm(self):
        self.assertTrue(self.qc.has_nterm_acetylation_ionquant("n[42.0106]AAAEEEDGGPEGPNR"))

    def test_lysine(self):
        self.assertFalse(self.qc.has_nterm_acetylation_ionquant("n[304.2071]PEPK[42.0106]R"))


if __name__ == "__main__":
    unittest.main()


class TestCombinePsmFiles(unittest.TestCase):
    """combine_psm_files replaces the module's head -1 / tail -n +2 concatenation of the groups' psm.tsv."""

    qc = load("tmt_qc")

    def test_header_once_rows_in_order(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            a, b = tmp / "a.tsv", tmp / "b.tsv"
            a.write_bytes(b"Spectrum\tPeptide\nA.1\tPEPK\n")
            b.write_bytes(b"Spectrum\tPeptide\nB.1\tKPEP\nB.2\tPEP")
            out = self.qc.combine_psm_files([a, b], tmp / "combined_psm.tsv")
            self.assertEqual(out.read_bytes(), b"Spectrum\tPeptide\nA.1\tPEPK\nB.1\tKPEP\nB.2\tPEP")

    def test_cli_takes_several_psm_files(self):
        parser_args = ["analyze", "--psm-file", "a.tsv", "b.tsv", "--tmt-type", "TMT16"]
        import argparse

        captured = {}
        original = argparse.ArgumentParser.parse_args

        def fake(self_parser, args=None, namespace=None):
            ns = original(self_parser, parser_args, namespace)
            captured["ns"] = ns
            raise SystemExit(0)

        argparse.ArgumentParser.parse_args = fake
        try:
            with self.assertRaises(SystemExit):
                self.qc.main()
        finally:
            argparse.ArgumentParser.parse_args = original
        self.assertEqual(captured["ns"].psm_file, ["a.tsv", "b.tsv"])
