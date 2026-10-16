"""End-to-end tests: real ``qpxc convert diann|fragpipe`` (qpx 1.1.4) + ``bin/qpx_export.py``.

Needs the qpx package (tests/python/requirements.txt).
Fixtures in tests/fixtures/qpx/ are tiny synthetic DIA-NN / FragPipe outputs with a
literal 0 and missing cells, a single-run DIA case, pipeline (``<source name>``)
and experiment annotation (``<experiment>_<bioreplicate>``) FragPipe naming, a fractionated FragPipe
case and a FragPipe table without MaxLFQ columns.

diann_multi: pg_matrix (PG.MaxLFQ) differs from the report PG.Quantity (2x),
Precursor.Normalised = 1.5 x Precursor.Quantity, P44444;P40000 is a two-member
group, P55555 is report-only (MaxLFQ 0) and P33333 x S2 fails Global.Q.Value
(0.03), so it is in qpx but not in pg_matrix.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pytest

NAN = np.nan

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "tests" / "fixtures" / "qpx"
EXPORT = REPO / "bin" / "qpx_export.py"
QPXC = Path(sys.executable).parent / "qpxc"
PREFIX = "T"
LINEAR = "NaN = not quantified (unstored cell, NA/empty, or source value <= 0; linear scale); dense float32"


def qpxc_convert(engine, case: Path, outdir: Path, pg_matrix=True, annotation=None):
    if engine == "diann":
        args = ["--report-path", case / "report.tsv", "--sdrf-file", case / "sdrf.tsv",
                "--diann-log", case / "diann.log.txt", "--qvalue-threshold", "0.05"]
        if pg_matrix:
            args += ["--pg-matrix-path", case / "pg_matrix.tsv"]
    else:
        args = ["--ion-file", case / "combined_ion.tsv", "--pg-file", case / "combined_protein.tsv",
                "--sdrf-file", case / "sdrf.tsv"] + (["--experiment-annotation-file", annotation] if annotation
                                                     else annotation_args(case))
    cmd = [QPXC, "convert", engine, *args, "--output-folder", outdir, "--output-prefix", PREFIX]
    result = subprocess.run(list(map(str, cmd)), capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return outdir


def annotation_args(case: Path):
    annotation = case / "experiment_annotation.tsv"
    return ["--experiment-annotation-file", annotation] if annotation.is_file() else []


def export(engine, case: Path, qpx_folder: Path, outdir: Path, sdrf=None, extra=(), pg_matrix=None, annotation=True):
    cmd = [sys.executable, EXPORT, "--engine", engine, "--sdrf-file", sdrf or case / "sdrf.tsv",
           "--qpx-folder", qpx_folder, "--prefix", PREFIX, "--outdir", outdir, *extra]
    if engine == "diann":
        cmd += ["--pg-matrix-file", pg_matrix or case / "pg_matrix.tsv"]
    if engine == "fragpipe":
        cmd += ["--pg-file", case / "combined_protein.tsv"] + (annotation_args(case) if annotation else [])
    return subprocess.run(list(map(str, cmd)), capture_output=True, text=True)


def last_error(result):
    return [line for line in result.stderr.splitlines() if line.startswith("ERROR:")][-1]


def read_summary(outdir):
    rows = dict(line.split("\t") for line in (outdir / f"{PREFIX}.summary.tsv").read_text().splitlines())
    stable = {k: v for k, v in rows.items() if not k.startswith("count.")}
    counts = {k: int(v) for k, v in rows.items() if k.startswith("count.")}
    return stable, counts


def stable_summary(engine, obs_key, sources):
    rows = {"engine": engine, "prefix": PREFIX, "obs_key": obs_key, "modalities": "precursors,proteins",
            "guards_passed": "G1,G2,G4,G5,G6,G7,G8,G9,G10,G16,G17" if engine == "diann"
            else "G1,G2,G3,G4,G5,G6,G7,G8,G9,G10,G11,G16"}
    for name, source in sources.items():
        if "." in name:
            rows[f"{name.split('.')[0]}.layers.raw.quantity_source"] = source
            continue
        rows[f"{name}.missing_value_convention"] = LINEAR
        rows[f"{name}.quantity_source"] = source
    return rows


DIANN_SOURCES = {
    "precursors": "Precursor.Normalised (DIA-NN report via qpx, pr_matrix filter Global.Q.Value & "
    "Global.PG.Q.Value <= 0.01; normalised)",
    "precursors.raw": "Precursor.Quantity (DIA-NN report via qpx; not normalised; pr_matrix cells only)",
    "proteins": "PG.MaxLFQ (DIA-NN pg_matrix; normalised)",
    "proteins.raw": "PG.Quantity (DIA-NN report via qpx; not normalised; pg_matrix cells only)",
}
FRAGPIPE_SOURCES = {"precursors": "<exp> Intensity (combined_ion.tsv)",
                    "proteins": "<exp> MaxLFQ Intensity (combined_protein.tsv)"}


@pytest.fixture(scope="module")
def converted(tmp_path_factory):
    """Run qpxc once per fixture case."""
    cache = {}

    def get(engine, name, **kwargs):
        key = (engine, name, tuple(sorted(kwargs.items())))
        if key not in cache:
            cache[key] = qpxc_convert(engine, DATA / name, tmp_path_factory.mktemp(name), **kwargs)
        return cache[key]

    return get


def assert_no_zero_and_nan_float32(path):
    adata = ad.read_h5ad(path)
    assert adata.X.dtype == np.float32
    assert isinstance(adata.X, np.ndarray), "X must be dense"
    assert not np.any(adata.X[~np.isnan(adata.X)] <= 0), "no 0 may remain in a linear matrix"
    return adata


class TestDiann:
    def test_multi_run_x_is_pg_matrix_and_raw_is_pg_quantity(self, converted, tmp_path):
        case = DATA / "diann_multi"
        result = export("diann", case, converted("diann", "diann_multi"), tmp_path)
        assert result.returncode == 0, result.stderr
        proteins = assert_no_zero_and_nan_float32(tmp_path / f"{PREFIX}.proteins.h5ad")
        assert list(proteins.obs_names) == ["S1", "S2", "S3"]
        # var = pg_matrix groups only: report-only P55555 (MaxLFQ 0) is not a var
        assert list(proteins.var_names) == ["P11111", "P22222", "P33333", "P44444"]
        expected = np.array([[1000, 2000, NAN, 4000], [1100, NAN, NAN, 4100], [1200, 2100, 3000, 4200]], np.float32)
        np.testing.assert_array_equal(proteins.X, expected)
        # raw = PG.Quantity (2x MaxLFQ) on the same cells; P33333 x S2 (Global.Q.Value 0.03) is masked
        np.testing.assert_array_equal(proteins.layers["raw"], expected * 2)
        assert proteins.uns["n_raw_cells_outside_pg_matrix"] == 1
        assert proteins.uns["n_report_only_groups_dropped"] == 1
        assert proteins.var.loc["P44444", "protein_group"] == "P44444;P40000"
        assert proteins.var.loc["P44444", "protein_accession"] == "P44444"
        assert not proteins.var["is_contaminant"].any()
        assert proteins.obs.loc["S2", "organism"] == "Homo sapiens"
        assert proteins.obs.loc["S2", "sdrf_biological_replicate"] == "2", "clashing SDRF column is prefixed"
        assert proteins.obs.loc["S2", "comment_file_uri"] == "s3://example-bucket/S2.raw"
        mdata = mu.read_h5mu(tmp_path / f"{PREFIX}.h5mu")
        assert list(mdata.mod) == ["precursors", "proteins"]
        precursors = mdata.mod["precursors"]
        assert list(precursors.obs_names) == list(proteins.obs_names)
        np.testing.assert_allclose(precursors.X, precursors.layers["raw"] * 1.5)  # Precursor.Normalised
        # pr_matrix filter: PEPTIDECR|3 x S2 (Global.Q.Value 0.03) is not a pr_matrix cell -> NaN in X and raw
        assert np.isnan(precursors.X[1, 2]) and np.isnan(precursors.layers["raw"][1, 2])
        assert int(np.isnan(precursors.X).sum()) == 6
        assert precursors.uns["n_raw_cells_outside_pr_matrix"] == 1
        assert "feature_mapping" in mdata.varp and mdata.varp["feature_mapping"].shape == (9, 9)
        stable, counts = read_summary(tmp_path)
        assert stable == stable_summary("diann", "run_stem", DIANN_SOURCES)
        assert counts["count.proteins.n_non_nan"] == 9 and counts["count.proteins.layers.raw.n_non_nan"] == 9

    def test_single_run(self, converted, tmp_path):
        result = export("diann", DATA / "diann_single", converted("diann", "diann_single"), tmp_path)
        assert result.returncode == 0, result.stderr
        proteins = assert_no_zero_and_nan_float32(tmp_path / f"{PREFIX}.proteins.h5ad")
        # pg_matrix: P11111, P22222, P44444; report-only P55555 (literal 0) is not a var
        assert proteins.shape == (1, 3)
        assert int(np.isnan(proteins.X).sum()) == 0
        assert proteins.obs.loc["S1", "sample_accession"] == "S1"
        stable, counts = read_summary(tmp_path)
        assert stable == stable_summary("diann", "run_stem", DIANN_SOURCES)
        assert counts["count.precursors.n_obs"] == 1

    def test_g2_without_pg_matrix(self, converted, tmp_path):
        folder = converted("diann", "diann_multi", pg_matrix=False)
        result = export("diann", DATA / "diann_multi", folder, tmp_path / "out")
        assert result.returncode == 3
        assert last_error(result) == (
            "ERROR: G2: no 'proteins' modality (built: ['precursors']); DIA-NN pg_matrix missing "
            "(--pg-matrix-path) or empty after the q-value filter"
        )
        assert not (tmp_path / "out").exists()

    def test_g17_pg_matrix_group_missing_from_qpx(self, converted, tmp_path):
        pg_matrix = tmp_path / "pg_matrix.tsv"
        text = (DATA / "diann_multi" / "pg_matrix.tsv").read_text()
        pg_matrix.write_text(text + "P77777\tP77777\tX_HUMAN\tX\tSynthetic\t10.0\t\t\n")
        result = export("diann", DATA / "diann_multi", converted("diann", "diann_multi"), tmp_path / "out",
                        pg_matrix=pg_matrix)
        assert result.returncode == 3
        assert last_error(result) == "ERROR: G17: pg_matrix protein groups ['P77777'] (1) are not in the qpx proteins view"
        assert not (tmp_path / "out").exists()

    def test_g16_duplicate_pg_rows_are_not_summed(self, converted, tmp_path):
        import duckdb

        folder = tmp_path / "qpx"
        shutil.copytree(converted("diann", "diann_multi"), folder)
        pg = folder / f"{PREFIX}.pg.parquet"
        with duckdb.connect() as con:
            con.execute(f"COPY (SELECT * FROM '{pg}' UNION ALL (SELECT * FROM '{pg}' LIMIT 1)) TO '{tmp_path}/pg.parquet'")
        shutil.move(tmp_path / "pg.parquet", pg)
        result = export("diann", DATA / "diann_multi", folder, tmp_path / "out")
        assert result.returncode == 3
        assert last_error(result) == (
            "ERROR: G16: pg.parquet has 13 rows for 12 distinct (run, protein) cells; qpx would sum the duplicates"
        )

    def test_g6_sdrf_from_other_experiment(self, converted, tmp_path):
        sdrf = tmp_path / "other.sdrf.tsv"
        sdrf.write_text((DATA / "diann_multi" / "sdrf.tsv").read_text().replace("S3.raw", "X9.raw"))
        result = export("diann", DATA / "diann_multi", converted("diann", "diann_multi"), tmp_path / "out", sdrf=sdrf)
        assert result.returncode == 3
        assert last_error(result) == (
            "ERROR: G6: 'precursors' obs_names differ from the expected set: missing=['X9'], unexpected=['S3']"
        )
        assert not (tmp_path / "out").exists()


class TestFragpipe:
    @pytest.mark.parametrize("case", ["fragpipe_pipeline", "fragpipe_annotation"])
    def test_precursors_relabelled_to_run_stems(self, converted, tmp_path, case):
        result = export("fragpipe", DATA / case, converted("fragpipe", case), tmp_path)
        assert result.returncode == 0, result.stderr
        proteins = assert_no_zero_and_nan_float32(tmp_path / f"{PREFIX}.proteins.h5ad")
        mdata = mu.read_h5mu(tmp_path / f"{PREFIX}.h5mu")
        precursors = mdata.mod["precursors"]
        runs = ["R1_DDA_S1", "R2_DDA_S2", "R3_DDA_S3"]
        assert list(proteins.obs_names) == runs and list(precursors.obs_names) == runs
        assert mdata.n_obs == 3, "modalities must share obs (not 6 unaligned rows)"
        assert list(precursors.obs["sample_accession"]) == ["DDA_S1", "DDA_S2", "DDA_S3"]
        # MaxLFQ 0 cells (DDA_S2 x GENEB, GENEC; DDA_S3 x GENEC) are NaN; all-zero contaminant is dropped
        assert int(np.isnan(proteins.X).sum()) == 3
        assert "contaminant_Q99999" not in proteins.var_names
        assert proteins.uns["n_source_zero_to_nan"] == 3, "source zeros are counted in combined_protein.tsv"
        assert proteins.uns["n_unstored_to_nan"] == 0
        np.testing.assert_array_equal(precursors.X[0], np.array([500, 900, 1500, 2000], dtype=np.float32))
        # contaminant_Q99999 (MaxLFQ 0 everywhere, Intensity 500) dropped; GENEC x S3 MaxLFQ 0 but Intensity 6100
        assert proteins.uns["n_proteins_all_missing_dropped"] == 1
        assert proteins.uns["n_dropped_with_intensity"] == 1
        assert proteins.uns["n_nan_with_intensity"] == 1
        assert proteins.var.loc["sp|P11111|GENEA_HUMAN", "protein_accession"] == "P11111"
        assert proteins.var.loc["sp|P11111|GENEA_HUMAN", "protein_group"] == "P11111;P11111-2"
        assert not proteins.var["is_contaminant"].any()
        assert proteins.obs.loc["R2_DDA_S2", "comment_data_file"] == "R2_DDA_S2.raw"
        if case == "fragpipe_annotation":
            assert proteins.obs.loc["R2_DDA_S2", "condition"] == "EXP"
            assert proteins.obs.loc["R2_DDA_S2", "replicate"] == "1"
        stable, counts = read_summary(tmp_path)
        assert stable == stable_summary("fragpipe", "run_stem", FRAGPIPE_SOURCES)
        assert counts["count.proteins.n_non_nan"] >= 9
        assert counts["count.proteins.n_nan_with_intensity"] == 1
        assert counts["count.proteins.n_proteins_all_missing_dropped"] == 1

    def test_fractionated_runs_collapse_to_experiments(self, converted, tmp_path):
        case = "fragpipe_fractionated"
        result = export("fragpipe", DATA / case, converted("fragpipe", case), tmp_path)
        assert result.returncode == 0, result.stderr
        proteins = assert_no_zero_and_nan_float32(tmp_path / f"{PREFIX}.proteins.h5ad")
        assert list(proteins.obs_names) == ["DDA_S1", "DDA_S2", "DDA_S3"]
        assert proteins.obs.loc["DDA_S1", "run_file_names"] == "R1_DDA_S1_F1;R1_DDA_S1_F2"
        assert int(np.isnan(proteins.X).sum()) == 3 and proteins.uns["n_source_zero_to_nan"] == 3
        mdata = mu.read_h5mu(tmp_path / f"{PREFIX}.h5mu")
        assert mdata.n_obs == 3
        assert list(mdata.mod["precursors"].obs_names) == list(proteins.obs_names)
        stable, _ = read_summary(tmp_path)
        assert stable == stable_summary("fragpipe", "experiment", FRAGPIPE_SOURCES)

    NO_MAXLFQ = (
        "ERROR: G3: combined_protein.tsv has no MaxLFQ Intensity / Total Intensity column for experiments "
        "['DDA_S1', 'DDA_S2', 'DDA_S3'] (only '<exp> Intensity'); qpx would silently drop the proteins "
        "modality. IonQuant writes MaxLFQ only with >= 3 experiment+bioreplicate combinations and "
        "ionquant.maxlfq=1; rerun with more experiments / ionquant.maxlfq=1, or disable the export with "
        "--skip_qpx_export"
    )

    def test_g3_precheck_without_maxlfq(self, tmp_path):
        case = DATA / "fragpipe_no_maxlfq"
        cmd = [sys.executable, EXPORT, "--engine", "fragpipe", "--precheck", "--pg-file",
               case / "combined_protein.tsv", "--sdrf-file", case / "sdrf.tsv"]
        result = subprocess.run(list(map(str, cmd)), capture_output=True, text=True)
        assert result.returncode == 3
        assert last_error(result) == self.NO_MAXLFQ

    def test_g3_export_without_maxlfq_even_after_qpxc_exit_0(self, converted, tmp_path):
        case = "fragpipe_no_maxlfq"
        result = export("fragpipe", DATA / case, converted("fragpipe", case), tmp_path / "out")
        assert result.returncode == 3
        assert last_error(result) == self.NO_MAXLFQ
        assert not (tmp_path / "out").exists()

    def test_precheck_passes_on_pipeline_table(self):
        case = DATA / "fragpipe_pipeline"
        cmd = [sys.executable, EXPORT, "--engine", "fragpipe", "--precheck", "--pg-file",
               case / "combined_protein.tsv", "--sdrf-file", case / "sdrf.tsv"]
        result = subprocess.run(list(map(str, cmd)), capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert "precheck OK" in result.stdout

    def test_g6_annotation_runs_not_in_sdrf(self, converted, tmp_path):
        case = tmp_path / "fragpipe_annotation"
        shutil.copytree(DATA / "fragpipe_annotation", case)
        (case / "sdrf.tsv").write_text((case / "sdrf.tsv").read_text().replace("R3_DDA_S3.raw", "R9_DDA_S3.raw"))
        result = export("fragpipe", case, converted("fragpipe", "fragpipe_annotation"), tmp_path / "out")
        assert result.returncode == 3
        assert last_error(result) == (
            "ERROR: G6: 'experiment->run mapping' obs_names differ from the expected set: missing=['R9_DDA_S3'], "
            "unexpected=['R3_DDA_S3']"
        )

    def test_g6_mzml_named_after_source_name(self, converted, tmp_path):
        """An annotation whose files are named after the source name, not the run file, fails G6."""
        case = DATA / "fragpipe_annotation"
        annotation = tmp_path / "experiment_annotation.tsv"
        text = (case / "experiment_annotation.tsv").read_text()
        for run, source in [("R1_DDA_S1", "DDA_S1"), ("R2_DDA_S2", "DDA_S2"), ("R3_DDA_S3", "DDA_S3")]:
            text = text.replace(f"/{run}.mzML", f"/{source}.mzML")
        annotation.write_text(text)
        qpx_folder = converted("fragpipe", "fragpipe_annotation", annotation=annotation)
        cmd = [sys.executable, EXPORT, "--engine", "fragpipe", "--sdrf-file", case / "sdrf.tsv",
               "--qpx-folder", qpx_folder, "--prefix", PREFIX, "--outdir", tmp_path / "out",
               "--pg-file", case / "combined_protein.tsv", "--experiment-annotation-file", annotation]
        result = subprocess.run(list(map(str, cmd)), capture_output=True, text=True)
        assert result.returncode == 3, result.stderr
        assert last_error(result) == (
            "ERROR: G6: 'experiment->run mapping' obs_names differ from the expected set: "
            "missing=['R1_DDA_S1', 'R2_DDA_S2', 'R3_DDA_S3'], unexpected=['DDA_S1', 'DDA_S2', 'DDA_S3']"
        )
        assert not (tmp_path / "out").exists()

    @pytest.mark.parametrize("case", ["fragpipe_pipeline", "fragpipe_annotation"])
    def test_proteins_h5ad_carries_dataset_level_uns(self, converted, tmp_path, case):
        result = export("fragpipe", DATA / case, converted("fragpipe", case), tmp_path)
        assert result.returncode == 0, result.stderr
        mdata = mu.read_h5mu(tmp_path / f"{PREFIX}.h5mu")
        proteins = ad.read_h5ad(tmp_path / f"{PREFIX}.proteins.h5ad")
        assert {"qpx_failed_modalities", "writer_version", "software_name"} <= set(mdata.uns)
        for key, value in mdata.uns.items():
            assert proteins.uns[key] == value, key



class TestCli:
    def test_help(self):
        result = subprocess.run([sys.executable, str(EXPORT), "--help"], capture_output=True, text=True)
        assert result.returncode == 0 and "--precheck" in result.stdout
