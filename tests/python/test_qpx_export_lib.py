"""Unit tests for ``bin/qpx_export_lib.py`` (no qpx needed).

Every guard has a failing case that asserts the exact ``ERROR:`` text and a
passing case. Missing values must be NaN, never 0
(docs/plans/2026-09-27_PLAN_qpx_h5ad_export.md, "Evidence for the missing-value rule").
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bin"))

import qpx_export_lib as lib
from qpx_export_lib import GuardError

DATA = Path(__file__).resolve().parents[1] / "fixtures" / "qpx"
NAN = np.nan


def raises(guard_call, message):
    """Assert that *guard_call* raises GuardError with exactly *message*."""
    with pytest.raises(GuardError) as excinfo:
        guard_call()
    assert str(excinfo.value) == message
    return excinfo.value


class TestRunFileStem:
    @pytest.mark.parametrize(
        "value, stem",
        [
            ("S1.raw", "S1"),
            ("/work/raw_files/R1.mzML", "R1"),
            ("run.d", "run"),
            ("x.mzML.gz", "x"),
            ("C:\\data\\A.RAW", "A"),
            ("no_suffix", "no_suffix"),
            ("run.d.zip", "run"),
            ("run.d.tar.gz", "run"),
        ],
    )
    def test_run_file_stem_matches_qpx_rule(self, value, stem):
        assert lib.run_file_stem(value) == stem

    @pytest.mark.parametrize(
        "value", ["S1.raw", "/w/R1.mzML", "x.mzML.gz", "run.d", "run.d.zip", "run.d.tar", "run.d.tgz", "a.wiff", "b.dia", "c.mgf", "no_suffix"]
    )
    def test_run_file_stem_is_the_installed_qpx_rule(self, value):
        from qpx.core.files import run_file_stem

        assert lib.run_file_stem(value) == run_file_stem(value)


class TestDensifyWithNan:
    def test_sparse_unstored_and_stored_zero_become_nan(self):
        # (0, 1) is an explicitly stored zero, as qpx stores DIA-NN PG.Quantity 0; (1, 0) is unstored
        matrix = sparse.csr_matrix(
            (np.array([1.0, 0.0, 2.0], dtype=np.float32), ([0, 0, 1], [0, 1, 1])), shape=(2, 2)
        )
        assert matrix.nnz == 3
        dense, stats = lib.densify_with_nan(matrix, linear=True)
        np.testing.assert_array_equal(dense, np.array([[1.0, NAN], [NAN, 2.0]], dtype=np.float32))
        assert dense.dtype == np.float32
        assert stats == {"n_unstored_to_nan": 1, "n_source_nan": 0, "n_source_zero_to_nan": 1, "n_nan": 2, "n_values": 4}

    def test_sparse_negative_value_becomes_nan_on_linear_scale(self):
        matrix = sparse.csr_matrix(np.array([[-3.0, 5.0]], dtype=np.float32))
        dense, stats = lib.densify_with_nan(matrix, linear=True)
        assert np.isnan(dense[0, 0]) and dense[0, 1] == 5.0
        assert stats["n_source_zero_to_nan"] == 1

    def test_log2_keeps_zero_and_only_nan_is_missing(self):
        dense, stats = lib.densify_with_nan(np.array([[0.0, NAN, -1.5]]), linear=False)
        np.testing.assert_array_equal(dense, np.array([[0.0, NAN, -1.5]], dtype=np.float32))
        assert stats == {"n_unstored_to_nan": 0, "n_source_nan": 1, "n_source_zero_to_nan": 0, "n_nan": 1, "n_values": 3}

    def test_g16_duplicate_sparse_cells_are_not_summed(self):
        # two stored entries for cell (0, 0): coo would sum them to 3.0
        matrix = sparse.coo_matrix((np.array([1.0, 2.0, 5.0]), ([0, 0, 1], [0, 0, 1])), shape=(2, 2))
        raises(
            lambda: lib.densify_with_nan(matrix, linear=True),
            "ERROR: G16: sparse matrix stores 1 duplicate (obs, var) cells; refusing to sum them",
        )

    def test_dense_input_is_not_modified(self):
        source = np.array([[0.0, 1.0]])
        lib.densify_with_nan(source, linear=True)
        assert source[0, 0] == 0.0

    def test_modality_uns_records_provenance(self):
        _, stats = lib.densify_with_nan(np.array([[0.0, 1.0]]), linear=True)
        uns = lib.modality_uns(stats, linear=True, quantity_source="PG.Quantity")
        assert uns == {
            "missing_value_convention": lib.LINEAR_CONVENTION,
            "n_source_zero_to_nan": 1,
            "n_unstored_to_nan": 0,
            "quantity_source": "PG.Quantity",
        }


class TestG1FailedModalities:
    def test_passes_when_record_is_empty(self):
        lib.check_failed_modalities({"qpx_failed_modalities": {}})

    def test_fails_when_record_is_missing(self):
        raises(
            lambda: lib.check_failed_modalities({}),
            "ERROR: G1: uns['qpx_failed_modalities'] is missing (qpx older than 1.1.3?); cannot verify the build",
        )

    def test_fails_when_a_modality_failed(self):
        raises(
            lambda: lib.check_failed_modalities({"qpx_failed_modalities": {"precursors": "ArrowInvalid: x"}}),
            "ERROR: G1: MuData modalities failed to build: {'precursors': 'ArrowInvalid: x'}",
        )


class TestG2G7ModalityPresent:
    def test_passes_when_present(self):
        lib.check_modality_present({"proteins": 1}, "proteins", "G2", "hint")

    def test_g2_fails_without_proteins(self):
        raises(
            lambda: lib.check_modality_present({"precursors": 1}, "proteins", "G2", "pg_matrix missing"),
            "ERROR: G2: no 'proteins' modality (built: ['precursors']); pg_matrix missing",
        )

    def test_g7_fails_without_precursors(self):
        raises(
            lambda: lib.check_modality_present({"proteins": 1}, "precursors", "G7", "no features"),
            "ERROR: G7: no 'precursors' modality (built: ['proteins']); no features",
        )


class TestG4NonEmpty:
    def test_passes_with_one_value(self):
        lib.check_nonempty("proteins", np.array([[NAN, 1.0]]))

    def test_fails_on_all_nan(self):
        raises(
            lambda: lib.check_nonempty("proteins", np.array([[NAN, NAN]])),
            "ERROR: G4: 'proteins' matrix is empty: shape (1, 2), 0 values",
        )

    def test_fails_on_zero_obs(self):
        raises(
            lambda: lib.check_nonempty("proteins", np.zeros((0, 3))),
            "ERROR: G4: 'proteins' matrix is empty: shape (0, 3), 0 values",
        )


class TestG5SampleAccession:
    def test_passes_when_filled(self):
        lib.check_sample_accession("proteins", pd.DataFrame({"sample_accession": ["S1"]}, index=["R1"]))

    def test_fails_without_column(self):
        raises(
            lambda: lib.check_sample_accession("proteins", pd.DataFrame(index=["R1"])),
            "ERROR: G5: 'proteins' obs has no sample_accession column: no SDRF sample matched any run",
        )

    def test_fails_on_empty_value(self):
        obs = pd.DataFrame({"sample_accession": ["S1", ""]}, index=["R1", "R2"])
        raises(
            lambda: lib.check_sample_accession("proteins", obs),
            "ERROR: G5: 'proteins' obs has empty sample_accession for ['R2']: runs not in the SDRF",
        )

    @pytest.mark.parametrize("empty", [None, pd.NA, np.nan, "None", "<NA>", "nan", "  "])
    def test_fails_on_none_na_nan_and_blank(self, empty):
        obs = pd.DataFrame({"sample_accession": pd.Series(["S1", empty], dtype=object, index=["R1", "R2"])})
        raises(
            lambda: lib.check_sample_accession("proteins", obs),
            "ERROR: G5: 'proteins' obs has empty sample_accession for ['R2']: runs not in the SDRF",
        )


class TestG6ObsNames:
    def test_passes_on_equal_sets(self):
        lib.check_obs_names("proteins", ["B", "A"], ["A", "B"])

    def test_fails_on_dropped_run(self):
        raises(
            lambda: lib.check_obs_names("proteins", ["A"], ["A", "B"]),
            "ERROR: G6: 'proteins' obs_names differ from the expected set: missing=['B'], unexpected=[]",
        )

    def test_fails_on_unexpected_run(self):
        raises(
            lambda: lib.check_obs_names("precursors", ["A", "Z"], ["A"]),
            "ERROR: G6: 'precursors' obs_names differ from the expected set: missing=[], unexpected=['Z']",
        )


class TestG8NanCount:
    def test_passes_on_match(self):
        lib.check_nan_count("proteins", np.array([[NAN, 1.0]]), 1)

    def test_fails_on_drift(self):
        raises(
            lambda: lib.check_nan_count("proteins", np.array([[NAN, 1.0]]), 0),
            "ERROR: G8: 'proteins' has 1 NaN but the source has 0 missing cells",
        )


class TestG9ValueRange:
    def test_passes_on_positive_linear(self):
        lib.check_value_range("proteins", np.array([[1.0, NAN]]), linear=True)

    def test_log2_accepts_zero(self):
        lib.check_value_range("protein", np.array([[0.0, -1.0]]), linear=False)

    def test_fails_on_zero_leak(self):
        raises(
            lambda: lib.check_value_range("proteins", np.array([[0.0, 1.0]]), linear=True),
            "ERROR: G9: 'proteins' contains 1 values <= 0 on a linear scale (zero leak)",
        )

    def test_fails_on_inf(self):
        raises(
            lambda: lib.check_value_range("protein", np.array([[np.inf, 1.0]]), linear=False),
            "ERROR: G9: 'protein' contains 1 infinite values",
        )


class TestG10AlignedObs:
    def test_passes_on_same_order(self):
        lib.check_aligned_obs({"precursors": pd.Index(["A", "B"]), "proteins": pd.Index(["A", "B"])})

    def test_fails_on_different_order(self):
        raises(
            lambda: lib.check_aligned_obs({"proteins": pd.Index(["A", "B"]), "precursors": pd.Index(["B", "A"])}),
            "ERROR: G10: 'precursors' obs_names ['B', 'A'] are not aligned with 'proteins' obs_names ['A', 'B']",
        )


class TestG11FractionRows:
    MAPPING = {"E1": ["E1_F1", "E1_F2"], "E2": ["E2_F1"]}
    NAMES = ["E1_F1", "E1_F2", "E2_F1"]

    def test_passes_on_identical_rows_with_nan(self):
        dense = np.array([[1.0, NAN], [1.0, NAN], [2.0, 3.0]], dtype=np.float32)
        lib.check_fraction_rows_identical(dense, self.NAMES, self.MAPPING)

    def test_fails_on_differing_fractions(self):
        dense = np.array([[1.0, NAN], [1.0, 5.0], [2.0, 3.0]], dtype=np.float32)
        raises(
            lambda: lib.check_fraction_rows_identical(dense, self.NAMES, self.MAPPING),
            "ERROR: G11: fraction rows of experiment 'E1' (['E1_F1', 'E1_F2']) differ; "
            "refusing to collapse them to one obs",
        )


class TestReadExperimentRuns:
    def test_pipeline_naming_uses_sdrf_source_name(self):
        mapping = lib.read_experiment_runs(DATA / "fragpipe_pipeline" / "sdrf.tsv")
        assert mapping == {"DDA_S1": ["R1_DDA_S1"], "DDA_S2": ["R2_DDA_S2"], "DDA_S3": ["R3_DDA_S3"]}

    def test_allinone_naming_uses_annotation(self):
        case = DATA / "fragpipe_allinone"
        mapping = lib.read_experiment_runs(case / "sdrf.tsv", case / "experiment_annotation.tsv")
        assert mapping == {"DDA_S1_1": ["R1_DDA_S1"], "DDA_S2_1": ["R2_DDA_S2"], "DDA_S3_1": ["R3_DDA_S3"]}

    def test_fractions_group_under_one_source_name(self):
        mapping = lib.read_experiment_runs(DATA / "fragpipe_fractionated" / "sdrf.tsv")
        assert mapping["DDA_S1"] == ["R1_DDA_S1_F1", "R1_DDA_S1_F2"]

    def test_fails_on_annotation_without_sample_column(self, tmp_path):
        bad = tmp_path / "experiment_annotation.tsv"
        bad.write_text("file\tsample_name\nR1.mzML\tS1\n")
        raises(
            lambda: lib.read_experiment_runs(DATA / "fragpipe_pipeline" / "sdrf.tsv", bad),
            "ERROR: G3: experiment annotation experiment_annotation.tsv lacks column(s) ['sample']",
        )

    def test_sdrf_run_stems_is_the_expected_run_set(self):
        assert lib.sdrf_run_stems(DATA / "diann_multi" / "sdrf.tsv") == ["S1", "S2", "S3"]


class TestG3FragpipePrecheck:
    def columns(self, case):
        return lib.read_tsv_header(DATA / case / "combined_protein.tsv")

    def test_passes_with_maxlfq_columns(self):
        mapping = lib.read_experiment_runs(DATA / "fragpipe_pipeline" / "sdrf.tsv")
        source = lib.precheck_fragpipe_columns(self.columns("fragpipe_pipeline"), mapping)
        assert source == "<exp> MaxLFQ Intensity (combined_protein.tsv)"

    def test_total_intensity_columns_are_the_primary_source(self):
        columns = ["A Spectral Count", "A Total Intensity", "A MaxLFQ Intensity"]
        source = lib.precheck_fragpipe_columns(columns, {"A": ["R1"]})
        assert source.startswith("<exp> Total Intensity")

    def test_fails_without_maxlfq_columns_and_names_the_cause(self):
        mapping = lib.read_experiment_runs(DATA / "fragpipe_no_maxlfq" / "sdrf.tsv")
        raises(
            lambda: lib.precheck_fragpipe_columns(self.columns("fragpipe_no_maxlfq"), mapping),
            "ERROR: G3: combined_protein.tsv has no MaxLFQ Intensity / Total Intensity column for experiments "
            "['DDA_S1', 'DDA_S2', 'DDA_S3'] (only '<exp> Intensity'); qpx would silently drop the proteins "
            "modality. IonQuant writes MaxLFQ only with >= 3 experiment+bioreplicate combinations and "
            "ionquant.maxlfq=1; rerun with more experiments / ionquant.maxlfq=1, or disable the export with "
            "--skip_qpx_export",
        )

    def test_fails_on_allinone_naming_without_annotation(self):
        mapping = lib.read_experiment_runs(DATA / "fragpipe_allinone" / "sdrf.tsv")
        raises(
            lambda: lib.precheck_fragpipe_columns(self.columns("fragpipe_allinone"), mapping),
            "ERROR: G3: experiments ['DDA_S1', 'DDA_S2', 'DDA_S3'] have no columns in combined_protein.tsv "
            "(experiments found in the header: ['DDA_S1_1', 'DDA_S2_1', 'DDA_S3_1']); allinone naming "
            "'<source name>_<bioreplicate>' needs --experiment-annotation-file",
        )

    def test_fails_on_header_experiment_missing_from_sdrf(self):
        raises(
            lambda: lib.precheck_fragpipe_columns(["A Spectral Count", "B Spectral Count"], {"A": ["R1"]}),
            "ERROR: G3: combined_protein.tsv experiments ['B'] are not in the SDRF/annotation",
        )


class TestFragpipeProteinSourceMissing:
    VARS = ["sp|P11111|GENEA_HUMAN", "sp|P22222|GENEB_HUMAN", "sp|P33333|GENEC_HUMAN", "sp|P44444|GENED_HUMAN"]
    EXPERIMENTS = ["DDA_S1", "DDA_S2", "DDA_S3"]

    def test_counts_source_zeros_of_kept_rows(self):
        counts = lib.fragpipe_protein_source_missing(
            DATA / "fragpipe_pipeline" / "combined_protein.tsv", self.EXPERIMENTS, self.VARS
        )
        # GENEB x S2, GENEC x S2, GENEC x S3; the all-zero contaminant row is not a var
        assert counts == {"n_source_zero_to_nan": 3, "n_unstored_to_nan": 0}

    def test_empty_cell_counts_as_absent(self, tmp_path):
        path = tmp_path / "combined_protein.tsv"
        path.write_text("Protein\tA MaxLFQ Intensity\tB MaxLFQ Intensity\nP1\t\t5.0\nP2\t0.0\t7.0\n")
        counts = lib.fragpipe_protein_source_missing(path, ["A", "B"], ["P1", "P2"])
        assert counts == {"n_source_zero_to_nan": 1, "n_unstored_to_nan": 1}

    def test_fails_when_rows_do_not_match_vars(self):
        raises(
            lambda: lib.fragpipe_protein_source_missing(
                DATA / "fragpipe_pipeline" / "combined_protein.tsv", self.EXPERIMENTS, self.VARS + ["sp|X|Y"]
            ),
            "ERROR: G8: combined_protein.tsv has 4 rows for the 5 protein vars (duplicate anchors?)",
        )


class TestCollapseAndAlign:
    MAPPING = {"E1": ["E1_F1", "E1_F2"], "E2": ["E2_F1", "E2_F2"]}

    def protein_obs(self):
        return pd.DataFrame(
            {"sample_accession": ["E1", "E1", "E2", "E2"], "fraction": ["1", "2", "1", "2"]},
            index=["E1_F1", "E1_F2", "E2_F1", "E2_F2"],
        )

    def test_collapse_keeps_shared_and_joins_differing_obs(self):
        dense = np.array([[1.0], [1.0], [NAN], [NAN]], dtype=np.float32)
        matrix, obs = lib.collapse_fragpipe_obs(dense, self.protein_obs(), self.MAPPING)
        np.testing.assert_array_equal(matrix, np.array([[1.0], [NAN]], dtype=np.float32))
        assert obs.to_dict("index") == {
            "E1": {"sample_accession": "E1", "fraction": "1;2", "run_file_names": "E1_F1;E1_F2"},
            "E2": {"sample_accession": "E2", "fraction": "1;2", "run_file_names": "E2_F1;E2_F2"},
        }

    def test_align_relabels_experiments_to_run_stems(self):
        protein_obs = pd.DataFrame({"sample_accession": ["S1", "S2"]}, index=["R1", "R2"])
        order, obs = lib.align_fragpipe_precursor_obs(["S2_1", "S1_1"], protein_obs,
                                                      {"S1_1": ["R1"], "S2_1": ["R2"]}, collapsed=False)
        assert order == [1, 0]
        assert list(obs.index) == ["R1", "R2"]

    def test_align_keeps_experiment_names_when_collapsed(self):
        protein_obs = pd.DataFrame({"sample_accession": ["E1", "E2"]}, index=["E1", "E2"])
        order, _ = lib.align_fragpipe_precursor_obs(["E2", "E1"], protein_obs, self.MAPPING, collapsed=True)
        assert order == [1, 0]


class TestSummary:
    def test_rows_are_sorted_and_counts_are_separate(self, tmp_path):
        dense, stats = lib.densify_with_nan(np.array([[0.0, 1.0]]), linear=True)
        rows = {"engine": "diann", **lib.summary_rows("", "proteins", dense, lib.modality_uns(stats, True, "X"))}
        path = tmp_path / "s.tsv"
        lib.write_summary_tsv(path, rows)
        lines = path.read_text().splitlines()
        assert lines == sorted(lines)
        assert "count.proteins.n_nan\t1" in lines
        assert "proteins.quantity_source\tX" in lines

    def test_rejects_tab_in_value(self, tmp_path):
        with pytest.raises(ValueError, match="contains a tab or newline"):
            lib.write_summary_tsv(tmp_path / "s.tsv", {"k": "a\tb"})



class TestObsMetadata:
    @pytest.mark.parametrize(
        "header, name",
        [
            ("characteristics[organism part]", "organism_part"),
            ("Characteristics[Organism]", "organism"),
            ("factor value[treatment dose]", "factor_treatment_dose"),
            ("comment[data file]", "comment_data_file"),
            ("comment[modification parameters].1", "comment_modification_parameters_1"),
            ("Source Name", "source_name"),
            ("assay name", "assay_name"),
        ],
    )
    def test_sanitize_column(self, header, name):
        assert lib.sanitize_column(header) == name

    def test_sdrf_table_is_indexed_by_run_stem(self):
        table = lib.sdrf_obs_table(DATA / "diann_multi" / "sdrf.tsv")
        assert list(table.index) == ["S1", "S2", "S3"]
        assert table.loc["S2", "biological_replicate"] == "2"
        assert table.loc["S1", "comment_modification_parameters_2"].startswith("NT=Acetyl")

    def test_sdrf_table_fails_on_run_listed_twice(self, tmp_path):
        lines = (DATA / "diann_multi" / "sdrf.tsv").read_text().splitlines()
        path = tmp_path / "dup.sdrf.tsv"
        path.write_text("\n".join(lines + [lines[1]]) + "\n")
        raises(lambda: lib.sdrf_obs_table(path), "ERROR: G5: SDRF dup.sdrf.tsv lists runs ['S1'] on more than one row")

    def test_sdrf_table_fails_on_sanitized_collision(self, tmp_path):
        path = tmp_path / "c.sdrf.tsv"
        path.write_text("comment[data file]\tcharacteristics[age]\tage\nS1.raw\t1\t2\n")
        raises(lambda: lib.sdrf_obs_table(path), "ERROR: G5: SDRF c.sdrf.tsv columns collide after sanitizing: ['age']")

    def test_join_prefixes_clashing_columns(self):
        obs = pd.DataFrame({"sample_accession": ["S1", "S2"], "biological_replicate": [1, 2]}, index=["R1", "R2"])
        extra = pd.DataFrame({"biological_replicate": ["1", "2"], "disease": ["x", "y"]}, index=["R2", "R1"])
        joined = lib.join_obs_metadata(obs, extra, "sdrf_")
        assert list(joined.columns) == ["sample_accession", "biological_replicate", "sdrf_biological_replicate", "disease"]
        assert joined.loc["R1", "disease"] == "y" and joined.loc["R1", "biological_replicate"] == 1

    def test_join_fails_on_run_without_metadata(self):
        obs = pd.DataFrame({"sample_accession": ["S1"]}, index=["R9"])
        raises(
            lambda: lib.join_obs_metadata(obs, pd.DataFrame({"x": ["1"]}, index=["R1"]), "sdrf_"),
            "ERROR: G5: obs ['R9'] have no metadata row (runs not in the SDRF/annotation)",
        )

    def test_annotation_table_has_condition_and_replicate(self):
        table = lib.annotation_obs_table(DATA / "fragpipe_allinone" / "experiment_annotation.tsv")
        assert table.to_dict("index")["R1_DDA_S1"] == {"condition": "EXP", "replicate": "1"}


class TestProteinIdentifiers:
    @pytest.mark.parametrize(
        "token, accession",
        [("sp|P12345|NAME_HUMAN", "P12345"), ("tr|A0A0|X_HUMAN", "A0A0"), ("P12345", "P12345"),
         ("contaminant_P00761", "contaminant_P00761"), ("contam_sp|P13645|K1C10_HUMAN", "P13645")],
    )
    def test_protein_accession(self, token, accession):
        assert lib.protein_accession(token) == accession

    @pytest.mark.parametrize(
        "token, expected",
        [("contaminant_P00761", True), ("contam_sp|P13645|K1C10_HUMAN", True), ("Cont_P1", True),
         ("CONTAM_P19001", True), ("sp|P12345|NAME_HUMAN", False), ("P12345", False)],
    )
    def test_is_contaminant(self, token, expected):
        assert lib.is_contaminant(token) is expected

    def test_protein_var_columns(self):
        table = lib.protein_var_columns(
            ["sp|P1|A_HUMAN", "P2"], [["sp|P1|A_HUMAN", "sp|P1-2|A_HUMAN", "sp|P1|A_HUMAN"], ["P2", "contaminant_P9"]]
        )
        assert table.to_dict("records") == [
            {"protein_accession": "P1", "protein_group": "P1;P1-2", "is_contaminant": False},
            {"protein_accession": "P2", "protein_group": "P2;contaminant_P9", "is_contaminant": True},
        ]


class TestPgMatrix:
    def test_read_indexes_by_canonical_membership_and_run_stem(self):
        values, groups = lib.read_pg_matrix(DATA / "diann_multi" / "pg_matrix.tsv")
        assert list(values.columns) == ["S1", "S2", "S3"]
        assert list(values.index) == ["P11111", "P22222", "P33333", "P40000;P44444"]
        assert groups["P40000;P44444"] == "P44444;P40000", "DIA-NN's own Protein.Group string is kept"
        assert np.isnan(values.loc["P22222", "S2"])

    def test_read_fails_on_repeated_group(self, tmp_path):
        path = tmp_path / "pg.tsv"
        path.write_text("Protein.Group\tProtein.Ids\tA.mzML\nP1;P2\tP1;P2\t1\nP2;P1\tP2;P1\t2\n")
        raises(lambda: lib.read_pg_matrix(path), "ERROR: G17: pg.tsv repeats protein groups ['P1;P2']")

    def matrix(self):
        return lib.read_pg_matrix(DATA / "diann_multi" / "pg_matrix.tsv")[0]

    def test_match_passes_on_identical_values(self):
        values = self.matrix()
        dense = values.to_numpy().T.astype(np.float32)
        lib.check_matches_pg_matrix(dense, ["S1", "S2", "S3"], list(values.index), values)

    def test_match_fails_on_extra_cell(self):
        values = self.matrix()
        dense = values.to_numpy().T.astype(np.float32)
        dense[1, 1] = 5.0  # S2 x P22222 is empty in pg_matrix
        raises(
            lambda: lib.check_matches_pg_matrix(dense, ["S1", "S2", "S3"], list(values.index), values),
            "ERROR: G17: proteins X has 10 values but pg_matrix has 9 (9 on the exported axes)",
        )

    def test_match_fails_on_dropped_group(self):
        values = self.matrix()
        dense = values.to_numpy().T.astype(np.float32)[:, :3]
        raises(
            lambda: lib.check_matches_pg_matrix(dense, ["S1", "S2", "S3"], list(values.index)[:3], values),
            "ERROR: G17: proteins X has 6 values but pg_matrix has 9 (6 on the exported axes)",
        )

    def test_match_fails_on_value_drift(self):
        values = self.matrix()
        dense = values.to_numpy().T.astype(np.float32)
        dense[0, 0] *= 1.001
        raises(
            lambda: lib.check_matches_pg_matrix(dense, ["S1", "S2", "S3"], list(values.index), values),
            "ERROR: G17: proteins X values differ from pg_matrix beyond float32 tolerance",
        )


class TestFragpipeDroppedProteins:
    VARS = TestFragpipeProteinSourceMissing.VARS

    def test_counts_dropped_rows_and_nan_cells_with_intensity(self):
        counts = lib.fragpipe_dropped_proteins(
            DATA / "fragpipe_pipeline" / "combined_protein.tsv", ["DDA_S1", "DDA_S2", "DDA_S3"], self.VARS
        )
        # contaminant_Q99999: MaxLFQ 0 everywhere, Intensity 500 in S1; GENEC x S3: MaxLFQ 0, Intensity 6100
        assert counts == {"n_proteins_all_missing_dropped": 1, "n_dropped_with_intensity": 1, "n_nan_with_intensity": 1}

    def test_fails_when_a_dropped_protein_had_maxlfq(self, tmp_path):
        path = tmp_path / "combined_protein.tsv"
        path.write_text("Protein\tA Intensity\tA MaxLFQ Intensity\nP1\t5\t4\nP2\t6\t3\n")
        raises(
            lambda: lib.fragpipe_dropped_proteins(path, ["A"], ["P1"]),
            "ERROR: G8: qpx dropped proteins with a positive MaxLFQ/Total intensity: ['P2']",
        )


class TestSummaryLayers:
    def test_layer_rows_and_integer_uns_counts(self):
        dense = np.array([[1.0, NAN]], dtype=np.float32)
        uns = {"missing_value_convention": "c", "quantity_source": "X", "n_extra": np.int64(4),
               "layers": {"raw": "R"}}
        rows = lib.summary_rows("", "proteins", dense, uns, {"raw": np.array([[NAN, NAN]])})
        assert rows["proteins.layers.raw.quantity_source"] == "R"
        assert rows["count.proteins.layers.raw.n_nan"] == "2"
        assert rows["count.proteins.n_extra"] == "4"
        assert rows["count.proteins.n_nan"] == "1"
