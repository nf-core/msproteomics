"""Unit tests for bin/philosopher_tables.py: synthetic pairs inside and outside the
measured Philosopher envelope (determinism doc, "Philosopher filter/report nondeterminism")."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

import golden_compare as gc  # noqa: E402
import philosopher_tables as pt  # noqa: E402

SIDE_E = gc.Side(["/runs/a"])
SIDE_A = gc.Side(["/runs/b"])


def psm_table(rows):
    head = "Spectrum\tSpectrum File\tPeptide\tCharge\tProbability\tQvalue\tProtein\tMapped Proteins"
    return "\n".join([head] + ["\t".join(map(str, r)) for r in rows]) + "\n"


def filter_log(threshold, n=100, decoy=1, parsing=("a", "b")):
    lines = [f'time="10:00:00" level=info msg="Parsing /runs/X/{p}.pep.xml"' for p in parsing]
    lines += [
        'time="10:00:00" level=info msg="2+ Charge profile" decoy=5 target=50',
        f'time="10:00:00" level=info msg="Converged to 1.00 % FDR with {n} PSMs" decoy={decoy} threshold={threshold} total={n + decoy}',
        'time="10:00:00" level=info msg="Converged to 0.95 % FDR with 9 Proteins" decoy=1 threshold=0.98 total=10',
        f'time="10:00:00" level=info msg="Final report numbers after FDR filtering, and post-processing" ions=5 peptides=4 proteins=3 psms={n}',
    ]
    return "\n".join(lines) + "\n"


# 10 PSMs; Probabilities strictly decreasing so T_above(row i) = i.
BASE = [
    [f"s{i}", "f.pep.xml", f"PEP{i % 3}", 2, f"{1 - i * 0.01:.4f}", f"{0.001 * i:.6f}", "sp|P1|A", "sp|P2|B"]
    for i in range(10)
]


class Base(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.e = Path(self._tmp.name) / "e"
        self.a = Path(self._tmp.name) / "a"
        self.e.mkdir()
        self.a.mkdir()

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, name, exp, act):
        (self.e / name).write_text(exp.replace("/runs/X", "/runs/a"))
        (self.a / name).write_text(act.replace("/runs/X", "/runs/b"))
        return self.e / name, self.a / name


class TestPsm(Base):
    def test_identical_matches(self):
        pe, pa = self.write("psm.tsv", psm_table(BASE), psm_table(BASE))
        self.assertIsNone(pt.compare_psm(pe, pa, SIDE_E, SIDE_A))

    def test_qvalue_within_bound(self):
        # row 4: T_above = 4, bound 1.01/4 = 0.2525
        rows = [r[:] for r in BASE]
        rows[4][5] = f"{float(rows[4][5]) + 0.25:.6f}"
        pe, pa = self.write("psm.tsv", psm_table(BASE), psm_table(rows))
        self.assertIsNone(pt.compare_psm(pe, pa, SIDE_E, SIDE_A))

    def test_qvalue_outside_bound(self):
        rows = [r[:] for r in BASE]
        rows[4][5] = f"{float(rows[4][5]) + 0.26:.6f}"
        pe, pa = self.write("psm.tsv", psm_table(BASE), psm_table(rows))
        self.assertIn("Qvalue", pt.compare_psm(pe, pa, SIDE_E, SIDE_A))

    def test_top_row_uses_max_t_above_1(self):
        rows = [r[:] for r in BASE]
        rows[0][5] = "1.0"  # T_above = 0 -> bound 1.01
        pe, pa = self.write("psm.tsv", psm_table(BASE), psm_table(rows))
        self.assertIsNone(pt.compare_psm(pe, pa, SIDE_E, SIDE_A))

    def test_other_column_must_match(self):
        rows = [r[:] for r in BASE]
        rows[3][2] = "OTHER"
        pe, pa = self.write("psm.tsv", psm_table(BASE), psm_table(rows))
        self.assertIn("Peptide", pt.compare_psm(pe, pa, SIDE_E, SIDE_A))

    def test_one_boundary_row_between_thresholds(self):
        extra = ["s10", "f.pep.xml", "PEPX", 2, "0.9050", "0.009", "sp|P1|A", ""]
        self.write("filter.log", filter_log(0.9050), filter_log(0.9052))
        pe, pa = self.write("psm.tsv", psm_table(BASE + [extra]), psm_table(BASE))
        self.assertIsNone(pt.compare_psm(pe, pa, SIDE_E, SIDE_A))

    def test_boundary_row_outside_thresholds(self):
        extra = ["s10", "f.pep.xml", "PEPX", 2, "0.8000", "0.009", "sp|P1|A", ""]
        self.write("filter.log", filter_log(0.9050), filter_log(0.9052))
        pe, pa = self.write("psm.tsv", psm_table(BASE + [extra]), psm_table(BASE))
        self.assertIn("not within thresholds", pt.compare_psm(pe, pa, SIDE_E, SIDE_A))

    def test_two_boundary_rows_rejected(self):
        extra = [["s10", "f", "PEPX", 2, "0.9051", "0.009", "p", ""], ["s11", "f", "PEPY", 2, "0.9051", "0.009", "p", ""]]
        self.write("filter.log", filter_log(0.9050), filter_log(0.9052))
        pe, pa = self.write("psm.tsv", psm_table(BASE + extra), psm_table(BASE))
        self.assertIn("max 1", pt.compare_psm(pe, pa, SIDE_E, SIDE_A))

    def test_boundary_row_without_filter_log_rejected(self):
        extra = ["s10", "f.pep.xml", "PEPX", 2, "0.9051", "0.009", "p", ""]
        pe, pa = self.write("psm.tsv", psm_table(BASE + [extra]), psm_table(BASE))
        self.assertIn("filter.log", pt.compare_psm(pe, pa, SIDE_E, SIDE_A))

    def test_row_order_matters(self):
        rows = [r[:] for r in BASE]
        rows[1], rows[2] = rows[2], rows[1]
        pe, pa = self.write("psm.tsv", psm_table(BASE), psm_table(rows))
        self.assertIn("order", pt.compare_psm(pe, pa, SIDE_E, SIDE_A))

    def test_duplicate_spectrum_keyed_with_file(self):
        rows = [r[:] for r in BASE] + [["s1", "g.pep.xml"] + BASE[1][2:]]
        pe, pa = self.write("psm.tsv", psm_table(rows), psm_table(rows))
        self.assertIsNone(pt.compare_psm(pe, pa, SIDE_E, SIDE_A))


def entry_table(kind, rows):
    heads = {
        "peptide.tsv": "Peptide\tPrev AA\tNext AA\tCharges\tProbability\tSpectral Count\tQvalue",
        "ion.tsv": "Peptide Sequence\tModified Sequence\tPrev AA\tNext AA\tCharge\tM/Z\tExpectation\tProbability\tSpectral Count\tQvalue",
        "protein.tsv": "Protein\tProtein ID\tTotal Spectral Count\tUnique Spectral Count\tRazor Spectral Count",
    }
    return "\n".join([heads[kind]] + ["\t".join(map(str, r)) for r in rows]) + "\n"


PEP = [[f"PEP{i}", "K", "A", "2", f"{1 - i * 0.1:.3f}", "3", f"{0.001 * i:.4f}"] for i in range(3)]
ION = [
    ["PEPTIDE", "", "K", "A", "2", "400.1", "0.01", "0.99", "2", "0.0"],
    ["MPEPTIDE", "M[147.0354]PEPTIDE", "K", "A", "2", "450.2", "0.02", "0.95", "1", "0.001"],
    ["LCIDIIK", "", "R", "S", "3", "300.3", "0.03", "0.90", "1", "0.002"],
]
PROT = [["sp|P1|A", "P1", "5", "4", "4"], ["sp|P2|B", "P2", "2", "2", "2"]]


class TestEntryTables(Base):
    def _with_boundary_psm(self, peptide="PEP1", charge=2, protein="sp|P1|A"):
        extra = ["s10", "f.pep.xml", peptide, charge, "0.9051", "0.009", protein, ""]
        self.write("filter.log", filter_log(0.9050), filter_log(0.9052))
        self.write("psm.tsv", psm_table(BASE + [extra]), psm_table(BASE))

    def test_spectral_count_moves_by_boundary_psm(self):
        self._with_boundary_psm()
        rows = [r[:] for r in PEP]
        rows[1][5] = "4"
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", PEP), entry_table("peptide.tsv", rows))
        self.assertIsNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))

    def test_ionquant_feature_columns_of_boundary_entry(self):
        # 12-run LFQ set, experiment_1_1: FragPipe's filter on identical inputs gives 15014 or 15016 PSMs (5/7 of 12
        # runs); IonQuant then quantifies the boundary PSM's peptide from a different PSM set
        head = "Peptide\tPrev AA\tNext AA\tCharges\tProbability\tSpectral Count\tIntensity\tQvalue"
        exp = [["PEP0", "K", "A", "2", "0.999", "3", "100.5", "0.0"], ["PEP1", "K", "A", "2", "0.95", "3", "25992.457", "0.001"]]
        act = [r[:] for r in exp]
        act[1][5], act[1][6] = "4", "38817.73"
        table = lambda rows: "\n".join([head] + ["\t".join(r) for r in rows]) + "\n"  # noqa: E731
        self._with_boundary_psm()
        pe, pa = self.write("peptide.tsv", table(exp), table(act))
        self.assertIsNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))
        act[0][6] = "101.0"  # no boundary PSM maps to PEP0
        pe, pa = self.write("peptide.tsv", table(exp), table(act))
        self.assertIn("Intensity", pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))

    def test_spectral_count_without_boundary_psm_rejected(self):
        rows = [r[:] for r in PEP]
        rows[1][5] = "4"
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", PEP), entry_table("peptide.tsv", rows))
        self.assertIn("boundary PSM", pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))

    def test_spectral_count_moved_by_two_rejected(self):
        self._with_boundary_psm()
        rows = [r[:] for r in PEP]
        rows[1][5] = "5"
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", PEP), entry_table("peptide.tsv", rows))
        self.assertIsNotNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))

    def test_protein_counts_follow_boundary_psm_protein(self):
        self._with_boundary_psm(protein="sp|P1|A")
        rows = [r[:] for r in PROT]
        rows[0][2] = "6"
        pe, pa = self.write("protein.tsv", entry_table("protein.tsv", PROT), entry_table("protein.tsv", rows))
        self.assertIsNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "protein.tsv", full=False))
        rows[1][2] = "3"
        pe, pa = self.write("protein.tsv", entry_table("protein.tsv", PROT), entry_table("protein.tsv", rows))
        self.assertIsNotNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "protein.tsv", full=False))

    def test_peptide_next_aa_and_qvalue_vary_in_dda_scope(self):
        # Measured: Philosopher filter+report on identical inputs gave Next AA "A" (5/6 runs) vs "L"
        # (1/6), and Qvalue 0.00123903926801 vs 0.00123974823574 (4/8 runs each)
        exp = [["LIEALLR", "K", "A", "2", "0.999", "3", "0.00123903926801"]] + [r[:] for r in PEP]
        act = [r[:] for r in exp]
        act[0][2] = "L"
        act[0][6] = "0.00123974823574"
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", exp), entry_table("peptide.tsv", act))
        self.assertIsNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))

    def test_peptide_qvalue_outside_bound_in_dda_scope(self):
        many = [[f"P{i}", "K", "A", "2", f"{1 - i * 0.01:.3f}", "3", "0.0010"] for i in range(30)]
        rows = [r[:] for r in many]
        rows[20][6] = "0.2"
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", many), entry_table("peptide.tsv", rows))
        self.assertIn("Qvalue", pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))

    def test_peptide_rows_without_probability(self):
        # IonQuant MBR rows in a rewritten peptide.tsv have empty Probability/Qvalue cells
        exp = [r[:] for r in PEP] + [["MBRPEP", "K", "A", "2", "", "0", ""]]
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", exp), entry_table("peptide.tsv", exp))
        self.assertIsNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))
        act = [r[:] for r in exp]
        act[1][6] = "0.0011"  # ranked row: inside the envelope
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", exp), entry_table("peptide.tsv", act))
        self.assertIsNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))
        act[3][6] = "0.01"  # unranked row: no envelope
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", exp), entry_table("peptide.tsv", act))
        self.assertIn("without a Philosopher Probability", pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))

    def test_peptide_other_columns_strict_in_dda_scope(self):
        rows = [r[:] for r in PEP]
        rows[1][3] = "3"  # Charges
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", PEP), entry_table("peptide.tsv", rows))
        self.assertIn("Charges", pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=False))

    def test_peptide_qvalue_rules_in_dia_scope(self):
        rows = [r[:] for r in PEP]
        rows[2][6] = "0.0"  # one side 0, other 0.002 <= 0.01
        rows[1][1] = "R"  # Prev AA dropped
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", PEP), entry_table("peptide.tsv", rows))
        self.assertIsNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=True))
        rows[1][6] = "0.9"  # T_above(row1)=1 -> bound 2.0, still inside
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", PEP), entry_table("peptide.tsv", rows))
        self.assertIsNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=True))

    def test_peptide_qvalue_outside_bound(self):
        many = [[f"P{i}", "K", "A", "2", f"{1 - i * 0.01:.3f}", "3", "0.0010"] for i in range(30)]
        rows = [r[:] for r in many]
        rows[20][6] = "0.2"  # bound 2/20 = 0.1, zero-rule does not apply
        pe, pa = self.write("peptide.tsv", entry_table("peptide.tsv", many), entry_table("peptide.tsv", rows))
        self.assertIn("Qvalue", pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "peptide.tsv", full=True))

    def test_ion_dia_rules(self):
        rows = [r[:] for r in ION]
        rows[0], rows[2] = rows[2], rows[0]  # sorted rows
        rows[1][6] = "429"  # Expectation dropped
        rows[1][2] = "R"  # Prev AA dropped
        pe, pa = self.write("ion.tsv", entry_table("ion.tsv", ION), entry_table("ion.tsv", rows))
        self.assertIsNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "ion.tsv", full=True))

    def test_ion_modified_sequence_as_mass_multiset(self):
        two = [["", "PEPM[147.0354]K[229.1629]", "K", "A", "2", "500.0", "0.1", "0.9", "1", "0.0"]]
        moved = [["", "PEPM[229.1629]K[147.0354]", "K", "A", "2", "500.0", "0.1", "0.9", "1", "0.0"]]
        other = [["", "PEPM[147.0354]K[147.0354]", "K", "A", "2", "500.0", "0.1", "0.9", "1", "0.0"]]
        pe, pa = self.write("ion.tsv", entry_table("ion.tsv", two), entry_table("ion.tsv", moved))
        self.assertIsNone(pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "ion.tsv", full=True))
        pe, pa = self.write("ion.tsv", entry_table("ion.tsv", two), entry_table("ion.tsv", other))
        self.assertIn("Modified Sequence", pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "ion.tsv", full=True))

    def test_ion_expectation_strict_in_dda_scope(self):
        rows = [r[:] for r in ION]
        rows[1][6] = "0.03"
        pe, pa = self.write("ion.tsv", entry_table("ion.tsv", ION), entry_table("ion.tsv", rows))
        self.assertIn("Expectation", pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "ion.tsv", full=False))


class TestFilterLog(Base):
    def test_parsing_order_and_counts_within_envelope(self):
        e = filter_log(0.7918, n=9803, decoy=98, parsing=("r2", "r4"))
        a = filter_log(0.7916, n=9805, decoy=99, parsing=("r4", "r2")).replace("psms=9805", "psms=9804")
        pe, pa = self.write("filter.log", e, a)
        self.assertIsNone(pt.compare_filter_log(pe, pa, SIDE_E, SIDE_A))

    def test_count_outside_envelope(self):
        pe, pa = self.write("filter.log", filter_log(0.79, n=100), filter_log(0.79, n=103))
        self.assertIn("count line", pt.compare_filter_log(pe, pa, SIDE_E, SIDE_A))

    def test_threshold_outside_envelope(self):
        pe, pa = self.write("filter.log", filter_log(0.7918), filter_log(0.7915))
        self.assertIsNotNone(pt.compare_filter_log(pe, pa, SIDE_E, SIDE_A))

    def test_decoy_outside_envelope(self):
        pe, pa = self.write("filter.log", filter_log(0.79, decoy=1), filter_log(0.79, decoy=3))
        self.assertIsNotNone(pt.compare_filter_log(pe, pa, SIDE_E, SIDE_A))

    def test_parsing_set_differs(self):
        pe, pa = self.write("filter.log", filter_log(0.79, parsing=("a", "b")), filter_log(0.79, parsing=("a", "c")))
        self.assertIn("Parsing", pt.compare_filter_log(pe, pa, SIDE_E, SIDE_A))

    def test_protein_and_charge_lines_exact(self):
        e = filter_log(0.79)
        pe, pa = self.write("filter.log", e, e.replace("decoy=5 target=50", "decoy=5 target=51"))
        self.assertIsNotNone(pt.compare_filter_log(pe, pa, SIDE_E, SIDE_A))
        pe, pa = self.write("filter.log", e, e.replace("with 9 Proteins", "with 10 Proteins").replace("total=10", "total=11"))
        self.assertIsNotNone(pt.compare_filter_log(pe, pa, SIDE_E, SIDE_A))

    def test_final_threshold_sequential_then_2d(self):
        seq = ['level=info msg="Converged to 1.00 % FDR with 10 PSMs" decoy=1 threshold=0.5 total=11',
               'level=info msg="9 PSMs" decoy=0 threshold=0.6 total=9']
        self.assertEqual(pt.final_psm_threshold(seq), 0.6)
        two = ['level=info msg="Converged to 1.00 % FDR with 10 PSMs" decoy=1 threshold=0.5 total=11',
               'level=info msg="Converged to 0.40 % FDR with 9 PSMs" decoy=0 threshold=0.7 total=9']
        self.assertEqual(pt.final_psm_threshold(two), 0.7)


if __name__ == "__main__":
    unittest.main()


class TestRunTables(Base):
    """Run-level tables: only entries a boundary PSM maps to may differ, only in its group's value columns
    (12-run LFQ set: FragPipe's filter on identical inputs gives 15014 or 15016 PSMs, 5 vs 7 of 12 runs)."""

    HEAD = "Peptide Sequence\tPrev AA\tG1 Spectral Count\tG1 Intensity\tG2 Spectral Count\tG2 Intensity"
    ROWS = [["PEP0", "K", "3", "100.5", "1", "7.0"], ["PEP1", "K", "2", "45829.277", "0", "0"]]

    def setUp(self):
        super().setUp()
        for side, thr, rows in ((self.e, 0.9050, BASE), (self.a, 0.9052, BASE + [["run1.180501.180501.2", "f.pep.xml", "PEP1", 2, "0.9051", "0.009", "sp|P1|A", ""]])):
            (side / "G1").mkdir()
            (side / "G1" / "filter.log").write_text(filter_log(thr))
            (side / "G1" / "psm.tsv").write_text(psm_table(rows))

    def table(self, rows, head=HEAD):
        return "\n".join([head] + ["\t".join(r) for r in rows]) + "\n"

    def compare(self, act_rows, name="combined_peptide.tsv", head=HEAD, exp_rows=ROWS):
        pe, pa = self.write(name, self.table(exp_rows, head), self.table(act_rows, head))
        return pt.compare_run_table(pe, pa, SIDE_E, SIDE_A, name)

    def test_boundary_entry_value_columns_may_differ(self):
        act = [r[:] for r in self.ROWS]
        act[1][2], act[1][3] = "3", "65584.39"
        self.assertIsNone(self.compare(act))

    def test_other_entry_must_be_exact(self):
        act = [r[:] for r in self.ROWS]
        act[1][3] = "65584.39"
        act[0][3] = "100.6"
        self.assertIn("no boundary PSM explains it", self.compare(act))

    def test_boundary_entry_other_columns_must_be_exact(self):
        act = [r[:] for r in self.ROWS]
        act[1][1] = "R"  # Prev AA
        self.assertIn("no boundary PSM explains it", self.compare(act))

    def test_boundary_entry_other_group_columns_must_be_exact(self):
        act = [r[:] for r in self.ROWS]
        act[1][5] = "1.0"  # G2 Intensity: the boundary PSM is in G1
        self.assertIn("no boundary PSM explains it", self.compare(act))

    def test_identical_tables_match_without_psm_siblings(self):
        (self.e / "G1" / "psm.tsv").unlink()
        self.assertIsNone(self.compare([r[:] for r in self.ROWS]))

    def test_msstats_only_the_boundary_psm_run(self):
        head = "ProteinName,PeptideSequence,PrecursorCharge,Run,Intensity"
        exp = ["sp|P1|A,PEP1,2,run1,45829.277", "sp|P1|A,PEP1,2,run2,NA", "sp|P1|A,PEP0,2,run1,10.0"]
        ok = [exp[0].replace("45829.277", "65584.39"), exp[1], exp[2]]
        bad = [exp[0], exp[1].replace("NA", "5.0"), exp[2]]
        for act, expect_none in ((ok, True), (bad, False)):
            pe, pa = self.write("msstats.csv", "\n".join([head] + exp) + "\n", "\n".join([head] + act) + "\n")
            r = pt.compare_run_table(pe, pa, SIDE_E, SIDE_A, "msstats.csv")
            self.assertEqual(r is None, expect_none, r)

    def test_reprint_protein_id_of_boundary_psm(self):
        head = "PROTID\tGENEID\tG1_SPC\tG2_SPC"
        exp = [["P1", "A", "6", "1"], ["P2", "B", "4", "4"]]
        act = [["P1", "A", "7", "1"], ["P2", "B", "4", "4"]]
        self.assertIsNone(self.compare(act, "reprint.spc.tsv", head, exp))
        act[1][2] = "5"  # P2 is a mapped protein only in expected rows; the boundary PSM maps to P1 only
        self.assertIsNotNone(self.compare(act, "reprint.spc.tsv", head, exp))


    def test_philosopher_trailing_tab_rows(self):
        # combined_protein.tsv as Philosopher abacus writes it (lib/aba/pro.go): no trailing tab on the
        # header, one on every data row; reprint.*.tsv: a short "na" row
        head = "Protein\tProtein ID\tG1 Spectral Count"
        exp = [["sp|P1|A", "P1", "6", ""], ["sp|P2|B", "P2", "4", ""]]
        self.assertIsNone(self.compare([r[:] for r in exp], "combined_protein.tsv", head, exp))
        act = [r[:] for r in exp]
        act[0][3] = "x"  # the extra field is compared too, even on the boundary PSM's protein P1
        self.assertIn("no boundary PSM explains it", self.compare(act, "combined_protein.tsv", head, exp))
        rep = [["na", "na", ""], ["P1", "A", "6", ""]]
        self.assertIsNone(self.compare([r[:] for r in rep], "reprint.spc.tsv", "PROTID\tGENEID\tG1_SPC\t", rep))

    def test_strict_table_rejects_ragged_rows(self):
        pe, _ = self.write("psm.tsv", "A\tB\n1\t2\t\n", "A\tB\n1\t2\n")
        with self.assertRaises(ValueError):
            pt.Table(pe)


class TestIonRepresentative(Base):
    """ion.tsv, every scope: Philosopher report on identical inputs (PXD001468 commonPTMs_DMO 6 runs,
    Basic-Search 4 runs) changes the row order and the representative PSM of an ion: Modified Sequence
    and Expectation are those of one of the ion's own PSMs."""

    PSM_HEAD = "Spectrum\tPeptide\tModified Peptide\tCharge\tExpectation\tProbability\tQvalue\tProtein"
    PSMS = [
        ["s1", "PEPMPK", "PEPM[147.0354]PK", "2", "0.001", "0.99", "0.0", "sp|P1|A"],
        ["s2", "PEPMPK", "PEPMP[113.0477]K", "2", "0.00002", "0.98", "0.0", "sp|P1|A"],
        ["s3", "AAAK", "", "2", "0.5", "0.97", "0.0", "sp|P1|A"],
    ]
    ION_HEAD = "Peptide Sequence\tModified Sequence\tPrev AA\tNext AA\tCharge\tM/Z\tExpectation\tProbability\tSpectral Count\tQvalue"

    def setUp(self):
        super().setUp()
        psm = "\n".join([self.PSM_HEAD] + ["\t".join(r) for r in self.PSMS]) + "\n"
        self.write("psm.tsv", psm, psm)

    def ions(self, rows):
        return "\n".join([self.ION_HEAD] + ["\t".join(r) for r in rows]) + "\n"

    EXP = [
        ["AAAK", "", "K", "A", "2", "200.1", "0.5", "0.97", "1", "0.0"],
        ["PEPMPK", "PEPM[147.0354]PK", "K", "A", "2", "400.2", "0.001", "0.99", "2", "0.0"],
        ["", "", "K", "A", "3", "300.3", "0.0000000361", "0.9", "0", "0.0"],
    ]

    def compare(self, act):
        pe, pa = self.write("ion.tsv", self.ions(self.EXP), self.ions(act))
        return pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, "ion.tsv", full=False)

    def test_other_representative_psm_and_row_order(self):
        act = [r[:] for r in self.EXP]
        act[1][1], act[1][6] = "PEPMP[113.0477]K", "0.00002"  # the ion's other PSM
        act[0], act[1] = act[1], act[0]
        self.assertIsNone(self.compare(act))

    def test_value_of_no_psm_of_the_ion_is_a_diff(self):
        act = [r[:] for r in self.EXP]
        act[1][6] = "0.5"  # AAAK's PSM, not PEPMPK's
        self.assertIn("not a value of the ion's own PSMs", self.compare(act))
        act = [r[:] for r in self.EXP]
        act[1][1] = "PEPM[147.0354]P[113.0477]K"
        self.assertIn("not a value of the ion's own PSMs", self.compare(act))

    def test_ion_without_psms_expectation_only(self):
        act = [r[:] for r in self.EXP]
        act[2][6] = "0.0000000000015"  # Spectral Count 0: measured to vary with no PSM to check against
        self.assertIsNone(self.compare(act))
        act[2][1] = "PEPK"
        self.assertIsNotNone(self.compare(act))

    def test_other_columns_stay_exact(self):
        act = [r[:] for r in self.EXP]
        act[1][7] = "0.95"  # Probability
        self.assertIn("Probability", self.compare(act))


def peptide_log(threshold, n=11720, decoy=118):
    return (
        f'time="10:00:00" level=info msg="Converged to 1.00 % FDR with {n} Peptides" decoy={decoy} '
        f"threshold={threshold} total={n + decoy}\n"
    )


class TestPeptideBoundary(Base):
    """Measured on PXD001468 commonPTMs_DMO: peptide threshold 0.772 vs 0.7725 moves one peptide (Probability 0.7720)."""

    EXTRA = ["GGSRGEEVGELSR", "K", "G", "2", "0.7720", "1", "0.0039"]

    def _compare(self, exp_rows, act_rows, log_e=None, log_a=None, kind="peptide.tsv"):
        if log_e is not None:
            self.write("filter.log", log_e, log_a)
        pe, pa = self.write(kind, entry_table(kind, exp_rows), entry_table(kind, act_rows))
        return pt.compare_entry_table(pe, pa, SIDE_E, SIDE_A, kind, full=False)

    def _logs(self, thr_e=0.772, thr_a=0.7725):
        return peptide_log(thr_e, n=11720, decoy=118), peptide_log(thr_a, n=11718, decoy=117)

    def test_one_peptide_on_lower_threshold_side_accepted(self):
        self.assertIsNone(self._compare(PEP + [self.EXTRA], PEP, *self._logs()))

    def test_one_peptide_on_lower_threshold_side_accepted_when_actual_is_lower(self):
        le, la = self._logs(0.7725, 0.772)
        self.assertIsNone(self._compare(PEP, PEP + [self.EXTRA], le, la))

    def test_extra_row_on_higher_threshold_side_rejected(self):
        self.assertIsNotNone(self._compare(PEP, PEP + [self.EXTRA], *self._logs()))

    def test_extra_row_with_equal_thresholds_rejected(self):
        self.assertIsNotNone(self._compare(PEP + [self.EXTRA], PEP, *self._logs(0.772, 0.772)))

    def test_peptide_outside_thresholds_rejected(self):
        extra = self.EXTRA[:4] + ["0.9000"] + self.EXTRA[5:]
        self.assertIsNotNone(self._compare(PEP + [extra], PEP, *self._logs()))

    def test_two_peptides_rejected(self):
        extra2 = ["PEPX", "K", "A", "2", "0.7722", "1", "0.0040"]
        self.assertIsNotNone(self._compare(PEP + [self.EXTRA, extra2], PEP, *self._logs()))

    def test_duplicate_key_extra_row_rejected(self):
        dup = PEP[0][:]
        self.assertIsNotNone(self._compare(PEP + [dup, self.EXTRA], PEP + [self.EXTRA], *self._logs()))
        self.assertIn("after dropping", self._compare(PEP + [dup], PEP, *self._logs()) or "")

    def test_row_count_difference_without_filter_log_rejected(self):
        self.assertIsNotNone(self._compare(PEP + [self.EXTRA], PEP))

    def test_ion_row_count_still_exact(self):
        self.assertIsNotNone(self._compare(ION, ION[:2], *self._logs(), kind="ion.tsv"))

    def test_protein_row_count_still_exact(self):
        self.assertIsNotNone(self._compare(PROT, PROT[:1], *self._logs(), kind="protein.tsv"))

    def test_peptide_log_threshold_tolerance_edge(self):
        for thr, ok in ((0.7725, True), (0.7726, False)):
            pe, pa = self.write("filter.log", peptide_log(0.772), peptide_log(thr, n=11718, decoy=117))
            res = pt.compare_filter_log(pe, pa, SIDE_E, SIDE_A)
            self.assertEqual(res is None, ok, f"{thr}: {res}")
            if not ok:
                self.assertIn("0.0005", res)

    def test_psm_log_threshold_still_tight(self):
        pe, pa = self.write("filter.log", filter_log(0.772), filter_log(0.7725))
        res = pt.compare_filter_log(pe, pa, SIDE_E, SIDE_A)
        self.assertIsNotNone(res)
        self.assertIn("0.0002", res)
