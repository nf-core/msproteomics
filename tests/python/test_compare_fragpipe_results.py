"""Unit tests for tests/tools/compare_fragpipe_results.py (stdlib unittest).

Run: make test-python
"""

import contextlib
import io
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import compare_fragpipe_results as cfr  # noqa: E402

N_ROWS = 20
EXPERIMENTS = ("exp_A", "exp_B")
TMT_CHANNELS = ("sample1", "sample2", "sample3")


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)] + ["\t".join(str(v) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n")


def peptide(i):
    return "PEPTIDE" + "ACDEFGHIK"[i % 9] * (1 + i // 9) + "K"


def build_fragpipe_dir(root, n_rows=N_ROWS):
    """Write a tiny FragPipe-style result directory with n_rows entities per table."""
    for e, exp in enumerate(EXPERIMENTS):
        d = root / exp
        psm_rows = []
        for i in range(n_rows):
            mods = "1M(15.9949), N-term(42.0106)" if i % 3 == 0 else ""
            psm_rows.append([f"{exp}.{i:05d}.{i:05d}.2", peptide(i), mods, f"{0.9 + i / 1000:.4f}",
                             10 + i, 1e6 * (i + 1) * (e + 1), "false", f"P{i % 7:05d}"])
        write_tsv(d / "psm.tsv", ["Spectrum", "Peptide", "Assigned Modifications", "Probability",
                                  "Hyperscore", "Intensity", "Is Decoy", "Protein ID"], psm_rows)
        write_tsv(d / "peptide.tsv", ["Peptide", "Intensity"],
                  [[peptide(i), 2e6 * (i + 1)] for i in range(n_rows)])
        write_tsv(d / "ion.tsv", ["Peptide Sequence", "Assigned Modifications", "Charge", "Intensity"],
                  [[peptide(i), "", 2, 3e6 * (i + 1)] for i in range(n_rows)])
        write_tsv(d / "protein.tsv", ["Protein ID", "Is Decoy", "Total Intensity"],
                  [[f"P{i:05d}", "false", 4e6 * (i + 1)] for i in range(n_rows)])
    quant_header = [f"{exp} {s}" for exp in EXPERIMENTS for s in ("Intensity", "MaxLFQ Intensity")]
    quant = [[5e5 * (i + 1) * (k + 1) for k in range(len(quant_header))] for i in range(n_rows)]
    write_tsv(root / "combined_protein.tsv", ["Protein ID"] + quant_header,
              [[f"P{i:05d}"] + q for i, q in enumerate(quant)])
    write_tsv(root / "combined_peptide.tsv", ["Peptide Sequence"] + quant_header,
              [[peptide(i)] + q for i, q in enumerate(quant)])
    write_tsv(root / "combined_modified_peptide.tsv", ["Peptide Sequence", "Assigned Modifications"] + quant_header,
              [[peptide(i), "N-term(42.0106)" if i % 2 else ""] + q for i, q in enumerate(quant)])
    write_tsv(root / "tmt-report" / "abundance_protein_MD.tsv",
              ["Index", "NumberPSM", "ReferenceIntensity", *TMT_CHANNELS],
              [[f"P{i:05d}", 3, 20.0, *(0.1 * i * (c + 1) for c in range(len(TMT_CHANNELS)))] for i in range(n_rows)])
    (root / "fragger.params").write_text(
        "database_name = /a/db.fasta # path\nnum_threads = 8\nprecursor_mass_lower = -20\n")


def run_cli(*argv):
    stderr = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(stderr):
        code = cfr.main([str(a) for a in argv])
    return code, stderr.getvalue()


class CompareFixture(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.ref = self.tmp / "fragpipe"
        self.test = self.tmp / "msproteomics"
        self.out = self.tmp / "out"
        build_fragpipe_dir(self.ref)
        build_fragpipe_dir(self.test)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def compare(self, *extra):
        return run_cli("--test", self.test, "--reference", self.ref, "--outdir", self.out, *extra)

    def report(self):
        return json.loads((self.out / "comparison.json").read_text())


class TestIdenticalDirectories(CompareFixture):
    def test_identical_dirs_pass_with_perfect_scores(self):
        code, _ = self.compare()
        self.assertEqual(code, cfr.EXIT_PASS)
        report = self.report()
        self.assertEqual(report["verdict"], "PASS")
        spearmans = [c["value"] for c in report["checks"] if c["check"].endswith("Spearman")]
        self.assertTrue(spearmans and all(v == 1.0 for v in spearmans), spearmans)
        self.assertEqual(report["combined"]["combined_protein.tsv"]["overlap"]["jaccard"], 1.0)

    def test_identical_dirs_write_markdown_report(self):
        self.compare()
        self.assertIn("**Verdict: PASS**", (self.out / "comparison.md").read_text())

    def test_tmt_channels_compared(self):
        self.compare()
        channels = self.report()["tmt"]["abundance_protein_MD.tsv"]["channels"]
        self.assertEqual(sorted(channels), sorted(TMT_CHANNELS))

    def test_nested_test_layout_matches_by_experiment_name(self):
        # msproteomics publishes experiments under philosopher/ and combined tables under ionquant/.
        shutil.rmtree(self.test)
        flat = self.tmp / "flat"
        build_fragpipe_dir(flat)
        for exp in EXPERIMENTS:
            shutil.move(str(flat / exp), str(self.test / "philosopher" / exp))
        for f in flat.glob("combined_*.tsv"):
            (self.test / "ionquant").mkdir(parents=True, exist_ok=True)
            shutil.move(str(f), str(self.test / "ionquant" / f.name))
        shutil.move(str(flat / "tmt-report"), str(self.test / "tmtintegrator"))
        code, err = self.compare()
        self.assertEqual(code, cfr.EXIT_PASS, err)


class TestDifferences(CompareFixture):
    def test_shifted_psm_counts_fail(self):
        shutil.rmtree(self.test / "exp_A")
        build_fragpipe_dir(self.tmp / "short", n_rows=N_ROWS - 2)
        shutil.move(str(self.tmp / "short" / "exp_A"), str(self.test / "exp_A"))
        code, _ = self.compare()
        self.assertEqual(code, cfr.EXIT_FAIL)
        failed = {c["check"] for c in self.report()["checks"] if c["status"] == "FAIL"}
        self.assertIn("exp_A: psms count deviation", failed)

    def test_relaxed_count_threshold_passes_counts(self):
        shutil.rmtree(self.test / "exp_A")
        build_fragpipe_dir(self.tmp / "short", n_rows=N_ROWS - 1)
        shutil.move(str(self.tmp / "short" / "exp_A"), str(self.test / "exp_A"))
        self.compare("--max-count-deviation", "0.1", "--min-jaccard", "0.9")
        failed = [c["check"] for c in self.report()["checks"] if c["status"] == "FAIL"]
        self.assertFalse([c for c in failed if "count deviation" in c], failed)

    def test_disjoint_spectrum_ids_fail(self):
        path = self.test / "exp_A" / "psm.tsv"
        path.write_text(path.read_text().replace("exp_A.", "x."))
        code, _ = self.compare()
        self.assertEqual(code, cfr.EXIT_FAIL)
        failed = {c["check"] for c in self.report()["checks"] if c["status"] == "FAIL"}
        self.assertIn("exp_A: psm:spectra Jaccard", failed)

    def test_scrambled_tmt_channel_fails_spearman(self):
        path = self.test / "tmt-report" / "abundance_protein_MD.tsv"
        lines = path.read_text().splitlines()
        rows = [ln.split("\t") for ln in lines[1:]]
        values = [r[3] for r in rows][::-1]
        for r, v in zip(rows, values):
            r[3] = v
        path.write_text("\n".join([lines[0]] + ["\t".join(r) for r in rows]) + "\n")
        code, _ = self.compare()
        self.assertEqual(code, cfr.EXIT_FAIL)
        failed = {c["check"] for c in self.report()["checks"] if c["status"] == "FAIL"}
        self.assertIn("abundance_protein_MD.tsv: channel sample1 Spearman", failed)

    def test_params_diff_reports_non_environment_keys(self):
        (self.test / "fragger.params").write_text(
            "database_name = /b/other.fasta\nnum_threads = 64\nprecursor_mass_lower = -10\n")
        self.compare()
        diff = self.report()["params"]["fragger.params"]
        self.assertEqual(list(diff["differing"]), ["precursor_mass_lower"])
        self.assertEqual(sorted(diff["differing_environment"]), ["database_name", "num_threads"])


class TestInputErrors(CompareFixture):
    def test_missing_required_column_errors_with_path(self):
        path = self.test / "exp_A" / "psm.tsv"
        path.write_text(path.read_text().replace("Assigned Modifications", "Mods"))
        code, err = self.compare()
        self.assertEqual(code, cfr.EXIT_ERROR)
        self.assertIn(str(path), err)
        self.assertIn("Assigned Modifications", err)

    def test_unmatched_experiment_names_error(self):
        for exp in EXPERIMENTS:
            (self.test / exp).rename(self.test / f"renamed_{exp}")
        code, err = self.compare()
        self.assertEqual(code, cfr.EXIT_ERROR)
        self.assertIn("--experiment-map", err)

    def test_experiment_map_resolves_renamed_experiments(self):
        for exp in EXPERIMENTS:
            (self.test / exp).rename(self.test / f"renamed_{exp}")
        code, err = self.compare("--experiment-map", *(f"renamed_{e}={e}" for e in EXPERIMENTS))
        self.assertEqual(code, cfr.EXIT_PASS, err)

    def test_combined_table_on_one_side_errors(self):
        (self.test / "combined_peptide.tsv").unlink()
        code, err = self.compare()
        self.assertEqual(code, cfr.EXIT_ERROR)
        self.assertIn("combined_peptide.tsv", err)

    def test_missing_directory_errors(self):
        code, err = run_cli("--test", self.tmp / "nope", "--reference", self.ref, "--outdir", self.out)
        self.assertEqual(code, cfr.EXIT_ERROR)
        self.assertIn("nope", err)


class TestHelpers(unittest.TestCase):
    def test_mod_key_is_order_and_precision_insensitive(self):
        a = cfr.mod_key("MPEPK", "1M(15.9949), N-term(42.0106)", 2, Path("x"))
        b = cfr.mod_key("MPEPK", "N-term(42.01060),1M(15.995)", 2, Path("x"))
        self.assertEqual(a, b)

    def test_mod_key_distinguishes_positions(self):
        self.assertNotEqual(cfr.mod_key("MPMK", "1M(15.9949)", 2, Path("x")),
                            cfr.mod_key("MPMK", "3M(15.9949)", 2, Path("x")))

    def test_mod_key_rejects_unparseable_token(self):
        with self.assertRaises(cfr.InputError):
            cfr.mod_key("PEPK", "garbage", 2, Path("x"))

    def test_spearman_handles_ties(self):
        self.assertAlmostEqual(cfr.spearman([1, 2, 2, 3], [10, 20, 20, 30]), 1.0)

    def test_quant_agreement_ignores_zeros_and_counts_one_sided(self):
        m = cfr.quant_agreement([(2.0, 4.0), (0.0, 1.0), (8.0, 8.0)], min_n=1, log=True)
        self.assertEqual((m["n"], m["only_reference"]), (2, 1))
        self.assertAlmostEqual(m["median_log2_ratio"], -0.5)


if __name__ == "__main__":
    unittest.main()
