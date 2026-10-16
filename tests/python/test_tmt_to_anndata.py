"""Tests for ``bin/tmt_to_anndata.py`` (TMT-Integrator abundance reports -> h5mu/h5ad).

Fixtures in tests/fixtures/qpx/tmt_* use the TMT-Integrator header layouts
(ReportInfo.java, ReportGenerator.writeTitle). The log2 set has NA, an empty
cell and a literal 0 (valid on log2 scale); the linear set has 0 (-> NaN).
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

import tmt_to_anndata as tmt

DATA = Path(__file__).resolve().parents[1] / "fixtures" / "qpx"
ANNOTATION = DATA / "tmt_annotation.tsv"
SAMPLES = ["P1_126", "P1_127N", "P1_127C", "P1_128N"]
LEVELS = {"gene": "gene", "protein": "protein", "peptide": "peptide",
          "single-site": "single_site", "multi-site": "multi_site"}


def run(tmp_path, files, log2="true", annotation=ANNOTATION):
    outdir = tmp_path / "out"
    code = tmt.main(["--abundance", *map(str, files), "--annotation", str(annotation),
                     "--log2transformed", log2, "--prefix", "T", "--outdir", str(outdir)])
    return code, outdir


def read_summary(outdir):
    rows = dict(line.split("\t") for line in (outdir / "T.tmt.summary.tsv").read_text().splitlines())
    return {k: v for k, v in rows.items() if not k.startswith("count.")}, {
        k: int(v) for k, v in rows.items() if k.startswith("count.")
    }


def expect_error(tmp_path, capsys, files, message, **kwargs):
    code, outdir = run(tmp_path, files, **kwargs)
    assert code == 3, "guard failures exit 3 (Nextflow terminates instead of retrying)"
    assert capsys.readouterr().err.strip() == message
    assert not outdir.exists(), "nothing may be written when a guard fails"


def copy_case(tmp_path, case):
    target = tmp_path / case
    shutil.copytree(DATA / case, target)
    return target


@pytest.fixture(scope="module")
def log2_out(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("log2")
    code, outdir = run(tmp, sorted((DATA / "tmt_log2").glob("abundance_*.tsv")))
    assert code == 0
    return outdir


class TestLog2Conversion:
    def test_writes_one_h5mu_per_norm_and_one_h5ad_per_level(self, log2_out):
        names = sorted(p.name for p in log2_out.iterdir())
        assert names == sorted(
            ["T.tmt.summary.tsv", "T.tmt_None.h5mu", "T.tmt_MD.h5mu", "T.tmt_gene_MD.h5ad", "T.tmt_protein_MD.h5ad"]
            + [f"T.tmt_{m}_None.h5ad" for m in LEVELS.values()]
        )

    def test_h5mu_has_one_modality_per_level(self, log2_out):
        mdata = mu.read_h5mu(log2_out / "T.tmt_None.h5mu")
        assert list(mdata.mod) == list(LEVELS.values())
        assert mdata.uns["prot_norm"] == "None" and bool(mdata.uns["log2transformed"]) is True

    @pytest.mark.parametrize("modality", list(LEVELS.values()))
    def test_h5ad_carries_dataset_level_uns(self, log2_out, modality):
        mdata = mu.read_h5mu(log2_out / "T.tmt_None.h5mu")
        adata = ad.read_h5ad(log2_out / f"T.tmt_{modality}_None.h5ad")
        assert set(mdata.uns) >= {"prot_norm", "source", "log2transformed_setting"}
        for key, value in mdata.uns.items():
            assert adata.uns[key] == value, key

    @pytest.mark.parametrize("level, modality", list(LEVELS.items()))
    def test_level_keeps_zero_and_turns_na_and_empty_into_nan(self, log2_out, level, modality):
        adata = ad.read_h5ad(log2_out / f"T.tmt_{modality}_None.h5ad")
        assert adata.X.dtype == np.float32
        assert list(adata.obs_names) == SAMPLES
        assert adata.n_vars == 3
        # row 1: 17.1 NA 16.9 0 ; row 2: 18.2 18.4 "" 18.0 ; row 3: 15.0 15.2 15.1 NA
        assert adata.X[3, 0] == 0.0, "0 is a valid log2 value and must be kept"
        assert np.isnan(adata.X[1, 0]) and np.isnan(adata.X[2, 1]) and np.isnan(adata.X[3, 2])
        assert int(np.isnan(adata.X).sum()) == 3
        assert adata.uns["groupby"] == level
        assert adata.uns["n_source_zero_to_nan"] == 0
        assert adata.uns["missing_value_convention"].startswith("NaN = not quantified (source NA or empty; log2")

    def test_protein_var_carries_annotation_and_refint_columns(self, log2_out):
        adata = ad.read_h5ad(log2_out / "T.tmt_protein_None.h5ad")
        assert list(adata.var_names) == ["P11111", "P21111", "P31111"]
        assert "Indistinguishable Proteins" in adata.var.columns and "RefInt_P1" in adata.var.columns
        assert "RefInt_P1" not in adata.obs_names

    def test_obs_joined_to_annotation(self, log2_out):
        adata = ad.read_h5ad(log2_out / "T.tmt_gene_MD.h5ad")
        assert adata.obs.loc["P1_127C", "channel"] == "127C"
        assert adata.obs.loc["P1_127C", "condition"] == "treat"

    def test_summary_rows(self, log2_out):
        stable, counts = read_summary(log2_out)
        assert stable["engine"] == "tmt"
        assert stable["norms"] == "MD,None"
        assert stable["guards_passed"] == "G4,G8,G9,G12,G13,G14,G15"
        assert stable["tmt_None.modalities"] == "gene,protein,peptide,single_site,multi_site"
        assert stable["tmt_MD.protein.quantity_source"] == "abundance_protein_MD.tsv"
        for modality in LEVELS.values():
            assert counts[f"count.tmt_None.{modality}.n_obs"] >= 4
            assert counts[f"count.tmt_None.{modality}.n_non_nan"] >= 1


class TestLinearConversion:
    def test_zero_and_na_become_nan(self, tmp_path):
        code, outdir = run(tmp_path, sorted((DATA / "tmt_linear").glob("abundance_*.tsv")), log2="false")
        assert code == 0
        adata = ad.read_h5ad(outdir / "T.tmt_protein_None.h5ad")
        # row 1: 120000.5 NA 98000 0 ; row 2: 5000 5100 "" 4900 ; row 3: 300 310 0.0 320
        assert int(np.isnan(adata.X).sum()) == 4
        assert np.isnan(adata.X[3, 0]) and np.isnan(adata.X[2, 2])
        assert not np.any(adata.X[~np.isnan(adata.X)] <= 0)
        assert adata.uns["n_source_zero_to_nan"] == 2
        _, counts = read_summary(outdir)
        assert counts["count.tmt_None.protein.n_source_zero_to_nan"] >= 2


class TestGuards:
    def test_g14_no_abundance_files(self):
        with pytest.raises(tmt.GuardError) as excinfo:
            tmt.classify_files([])
        assert str(excinfo.value) == "ERROR: G14: TMT-Integrator produced no abundance_*.tsv files"

    def test_g14_unglobbed_pattern(self, tmp_path, capsys):
        expect_error(
            tmp_path, capsys, ["abundance_*.tsv"],
            "ERROR: G14: 'abundance_*.tsv' is not a readable abundance_<level>_<norm>.tsv with level in "
            "['gene', 'multi-site', 'peptide', 'protein', 'single-site'] and norm in None/MD/GN/SL+IRS",
        )

    def test_g14_unknown_level(self, tmp_path, capsys):
        path = tmp_path / "abundance_modified-peptide_None.tsv"
        shutil.copy(DATA / "tmt_log2" / "abundance_peptide_None.tsv", path)
        expect_error(
            tmp_path, capsys, [path],
            f"ERROR: G14: '{path}' is not a readable abundance_<level>_<norm>.tsv with level in "
            "['gene', 'multi-site', 'peptide', 'protein', 'single-site'] and norm in None/MD/GN/SL+IRS",
        )

    def test_g12_wrong_header_layout(self, tmp_path, capsys):
        case = copy_case(tmp_path, "tmt_log2")
        path = case / "abundance_gene_None.tsv"
        text = path.read_text().replace("NumberPSM\t", "", 1)
        path.write_text(text)
        expect_error(
            tmp_path, capsys, [path],
            "ERROR: G12: abundance_gene_None.tsv: header ['Index', 'ProteinID', 'MaxPepProb', 'ReferenceIntensity', "
            "'P1_126'] does not match the TMT-Integrator gene layout "
            "['Index', 'NumberPSM', 'ProteinID', 'MaxPepProb', 'ReferenceIntensity']",
        )

    def test_g12_sample_column_not_in_annotation(self, tmp_path, capsys):
        annotation = tmp_path / "annotation.tsv"
        annotation.write_text(ANNOTATION.read_text().replace("P1_128N\tP1_128N", "P1_128N\tOTHER"))
        expect_error(
            tmp_path, capsys, [DATA / "tmt_log2" / "abundance_gene_None.tsv"],
            "ERROR: G12: abundance_gene_None.tsv: sample columns ['P1_128N'] are not annotation sample_name "
            "values ['Bridge', 'OTHER', 'P1_126', 'P1_127C', 'P1_127N']",
            annotation=annotation,
        )

    def test_g12_annotation_without_sample_name(self, tmp_path, capsys):
        annotation = tmp_path / "annotation.tsv"
        annotation.write_text("channel\tsample\n126\tA\n")
        expect_error(
            tmp_path, capsys, [DATA / "tmt_log2" / "abundance_gene_None.tsv"],
            "ERROR: G12: annotation annotation.tsv lacks column(s) ['plex', 'sample_name', 'condition', 'replicate']",
            annotation=annotation,
        )

    def test_g12_duplicated_index(self, tmp_path, capsys):
        case = copy_case(tmp_path, "tmt_log2")
        path = case / "abundance_gene_None.tsv"
        lines = path.read_text().splitlines()
        path.write_text("\n".join(lines + [lines[1]]) + "\n")
        expect_error(tmp_path, capsys, [path],
                     "ERROR: G12: abundance_gene_None.tsv: duplicated Index values ['GENE1']")

    def test_g12_non_numeric_cell(self, tmp_path, capsys):
        case = copy_case(tmp_path, "tmt_log2")
        path = case / "abundance_gene_None.tsv"
        path.write_text(path.read_text().replace("\t16.9\t", "\tabc\t", 1))
        expect_error(tmp_path, capsys, [path],
                     "ERROR: G12: abundance_gene_None.tsv: 1 sample cells are neither numeric nor NA/empty")

    def test_g13_samples_differ_between_levels(self, tmp_path, capsys):
        case = copy_case(tmp_path, "tmt_log2")
        path = case / "abundance_gene_None.tsv"
        lines = path.read_text().splitlines()
        # drop the last sample column (P1_128N) from the gene level only
        path.write_text("\n".join("\t".join(line.split("\t")[:-1]) for line in lines) + "\n")
        expect_error(
            tmp_path, capsys, [path, case / "abundance_protein_None.tsv"],
            "ERROR: G13: norm None: 'protein' samples ['P1_126', 'P1_127N', 'P1_127C', 'P1_128N'] differ from "
            "'gene' samples ['P1_126', 'P1_127N', 'P1_127C']",
        )

    def test_g15_log2_declared_on_linear_report(self, tmp_path, capsys):
        expect_error(
            tmp_path, capsys, [DATA / "tmt_linear" / "abundance_gene_None.tsv"],
            "ERROR: G15: abundance_gene_None.tsv: sample columns expected log2 (None: follows "
            "log2transformed=true (ReportGenerator.writeRatioAbnValues)) but values reach 120000; the report is linear",
        )

    def test_g9_infinite_value(self, tmp_path, capsys):
        case = copy_case(tmp_path, "tmt_log2")
        path = case / "abundance_gene_None.tsv"
        path.write_text(path.read_text().replace("\t16.9\t", "\tinf\t", 1))
        expect_error(tmp_path, capsys, [path], "ERROR: G9: 'tmt_None.gene' contains 1 infinite values")

    def test_g4_all_missing_report(self, tmp_path, capsys):
        case = copy_case(tmp_path, "tmt_log2")
        path = case / "abundance_gene_None.tsv"
        lines = path.read_text().splitlines()
        blank = [lines[0]] + ["\t".join(line.split("\t")[:5] + ["NA"] * 4) for line in lines[1:]]
        path.write_text("\n".join(blank) + "\n")
        expect_error(tmp_path, capsys, [path], "ERROR: G4: 'gene' matrix is empty: shape (4, 3), 0 values")


class TestScaleRule:
    """Per-norm scale from TMT-Integrator (PsmNormalizer.java:120-141, ReportGenerator.java:316-362)."""

    def test_sl_irs_is_log2_even_with_log2transformed_false(self, tmp_path):
        code, outdir = run(tmp_path, [DATA / "tmt_slirs" / "abundance_protein_SL+IRS.tsv"], log2="false")
        assert code == 0
        adata = ad.read_h5ad(outdir / "T.tmt_protein_SL+IRS.h5ad")
        # row 1: 17.3 NA 16.8 0.0 ; row 3: 8.2 8.3 -0.2 8.4 -> 0 and negative values are kept on log2
        assert adata.X[3, 0] == 0.0 and adata.X[2, 2] == np.float32(-0.2)
        assert int(np.isnan(adata.X).sum()) == 2
        assert adata.uns["scale"] == "log2" and bool(adata.uns["log2transformed"]) is True
        assert bool(adata.uns["log2transformed_setting"]) is False
        assert adata.uns["reference_intensity_scale"] == "linear"
        stable, _ = read_summary(outdir)
        assert stable["tmt_SL+IRS.protein.scale"] == "log2"
        assert stable["tmt_SL+IRS.protein.scale_rule"].startswith("SL+IRS: log2-like regardless of log2transformed")
        assert bool(mu.read_h5mu(outdir / "T.tmt_SL+IRS.h5mu").uns["log2transformed"]) is True

    def test_linear_declared_on_log2_report_fails(self, tmp_path, capsys):
        expect_error(
            tmp_path, capsys, [DATA / "tmt_log2" / "abundance_gene_MD.tsv"],
            "ERROR: G15: abundance_gene_MD.tsv: sample columns expected linear (MD: follows log2transformed=false "
            "(ReportGenerator.writeRatioAbnValues)) but the largest value is 18.3; the report is log2",
            log2="false",
        )

    def test_sl_irs_with_linear_values_fails(self, tmp_path, capsys):
        path = tmp_path / "abundance_protein_SL+IRS.tsv"
        shutil.copy(DATA / "tmt_linear" / "abundance_protein_None.tsv", path)
        expect_error(
            tmp_path, capsys, [path],
            "ERROR: G15: abundance_protein_SL+IRS.tsv: sample columns expected log2 (SL+IRS: log2-like regardless of "
            "log2transformed (PsmNormalizer.java:120-141, ReportGenerator.java:344-362)) but values reach 120000; "
            "the report is linear",
            log2="false",
        )

    def test_log2_reference_intensity_in_linear_report_fails(self, tmp_path, capsys):
        path = tmp_path / "abundance_gene_None.tsv"
        path.write_text((DATA / "tmt_linear" / "abundance_gene_None.tsv").read_text().replace("185363.8", "17.5"))
        expect_error(
            tmp_path, capsys, [path],
            "ERROR: G15: abundance_gene_None.tsv: ReferenceIntensity expected linear (None: follows "
            "log2transformed=false (ReportGenerator.writeRatioAbnValues)) but the largest value is 17.5; "
            "the report is log2",
            log2="false",
        )


class TestVarColumns:
    def test_numeric_header_columns_are_numbers(self, log2_out):
        var = ad.read_h5ad(log2_out / "T.tmt_multi_site_None.h5ad").var
        for column in ("Start", "End", "MaxPepProb", "Spectrum Number", "ReferenceIntensity"):
            assert var[column].dtype == np.float64, column
        protein = ad.read_h5ad(log2_out / "T.tmt_protein_None.h5ad").var
        assert protein["NumberPSM"].dtype == np.float64 and protein["RefInt_P1"].dtype == np.float64
        assert protein.loc["P11111", "NumberPSM"] == 3.0

    def test_protein_identifier_columns(self, log2_out):
        var = ad.read_h5ad(log2_out / "T.tmt_protein_None.h5ad").var
        assert var.loc["P11111", "protein_accession"] == "P11111"
        assert var.loc["P11111", "protein_group"] == "P11111;P11112;contaminant_P00761"
        assert bool(var.loc["P11111", "is_contaminant"]) is True
        assert bool(var.loc["P21111", "is_contaminant"]) is False
        gene = ad.read_h5ad(log2_out / "T.tmt_gene_None.h5ad").var
        assert gene.loc["GENE1", "protein_accession"] == "P11111"

    def test_g12_non_numeric_var_column(self, tmp_path, capsys):
        case = copy_case(tmp_path, "tmt_log2")
        path = case / "abundance_gene_None.tsv"
        path.write_text(path.read_text().replace("\t0.99\t", "\thigh\t", 1))
        expect_error(tmp_path, capsys, [path],
                     "ERROR: G12: abundance_gene_None.tsv: column 'MaxPepProb' has 1 values that are neither "
                     "numeric nor NA/empty")
