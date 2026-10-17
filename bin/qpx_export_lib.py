#!/usr/bin/env python3
"""Shared guards and matrix helpers for the qpx h5mu/h5ad export.

Pure library: depends on numpy, pandas, scipy and anndata only (no qpx import),
so it is unit-tested without qpx installed. Used by ``qpx_export.py`` (DIA-NN and
FragPipe via qpx) and ``tmt_to_anndata.py`` (TMT-Integrator reports).

Every guard raises :class:`GuardError` with a message ``ERROR: G<n>: ...``.
The CLIs turn it into exit code :data:`GUARD_EXIT_CODE` (3), distinct from a
Python crash (1) and an argparse error (2), so Nextflow can terminate instead of
retrying a deterministic guard failure. Nothing is written before all guards pass.

Missing-value rule: missing is
NaN, never 0. On a linear scale an unstored cell or a value <= 0 becomes NaN; on
a log2 scale only NA/empty becomes NaN (0 is a valid log2 value).
"""

from __future__ import annotations

import csv
import numbers
import re
from pathlib import Path, PurePosixPath

import numpy as np
import pandas as pd
from scipy import sparse

LINEAR_CONVENTION = (
    "NaN = not quantified (unstored cell, NA/empty, or source value <= 0; linear scale); dense float32"
)
LOG2_CONVENTION = "NaN = not quantified (source NA or empty; log2 scale, 0 is a valid value); dense float32"

# Same suffix rule as qpx.core.files.run_file_stem (qpx 1.1.5), so run names match qpx obs_names.
_RUN_FILE_SUFFIX = re.compile(r"(?i)\.(?:mzml(?:\.gz)?|mzxml|raw|d(?:\.(?:zip|tar(?:\.gz)?|tgz))?|wiff|mgf|dia)$")

GUARD_EXIT_CODE = 3

# Contaminant accession prefixes, compared case-insensitively:
# - "contaminant_": our FASTA contaminant entries (data/contaminants.fasta, conf/dia_base.config
#   contaminant_pattern), which DIA-NN and FragPipe report unchanged;
# - "contam_" / "Cont_": Philosopher / FragPipe (data/philosopher/lib/rep/protein.go:268,
#   data/FragPipe/.../Version.java:494 "Add `contam_` prefix to contaminants").
CONTAMINANT_PREFIXES = ("contaminant_", "contam_", "cont_")
# Tokens that mean "no value" in obs metadata (G5) after str(): None, pd.NA, NaN, empty.
EMPTY_TOKENS = {"", "none", "<na>", "nan"}

FRAGPIPE_MAXLFQ_HINT = (
    "IonQuant writes MaxLFQ only with >= 3 experiment+bioreplicate combinations and ionquant.maxlfq=1; "
    "rerun with more experiments / ionquant.maxlfq=1, or disable the export with --skip_qpx_export"
)


class GuardError(Exception):
    """A silent-failure guard fired. The message is the full ``ERROR: ...`` line."""

    def __init__(self, guard: str, message: str):
        self.guard = guard
        super().__init__(f"ERROR: {guard}: {message}")


def run_file_stem(value: str) -> str:
    """Return a run basename without its acquisition-file suffix (qpx rule)."""
    name = PurePosixPath(str(value).strip().replace("\\", "/")).name
    return _RUN_FILE_SUFFIX.sub("", name)


# ---------------------------------------------------------------------------
# Matrix conversion
# ---------------------------------------------------------------------------


def densify_with_nan(matrix, linear: bool) -> tuple[np.ndarray, dict[str, int]]:
    """Return a dense float32 copy of *matrix* with every missing cell as NaN.

    Unstored sparse cells become NaN. On a linear scale stored values <= 0 also
    become NaN; on a log2 scale they are kept. Stats count what was converted.
    """
    if sparse.issparse(matrix):
        coo = sparse.coo_matrix(matrix)
        n_cells = np.unique(coo.row.astype(np.int64) * coo.shape[1] + coo.col).size
        if n_cells != coo.nnz:
            raise GuardError("G16", f"sparse matrix stores {coo.nnz - n_cells} duplicate (obs, var) cells; refusing to sum them")
        dense = np.full(coo.shape, np.nan, dtype=np.float32)
        dense[coo.row, coo.col] = coo.data
        n_unstored = coo.shape[0] * coo.shape[1] - coo.nnz
    else:
        dense = np.array(matrix, dtype=np.float32, copy=True)
        n_unstored = 0
    n_source_nan = int(np.isnan(dense).sum())
    nonpositive = (dense <= 0) if linear else np.zeros(dense.shape, dtype=bool)
    dense[nonpositive] = np.nan
    stats = {
        "n_unstored_to_nan": int(n_unstored),
        "n_source_nan": n_source_nan - int(n_unstored),
        "n_source_zero_to_nan": int(nonpositive.sum()),
        "n_nan": int(np.isnan(dense).sum()),
        "n_values": int(dense.size),
    }
    return dense, stats


def modality_uns(stats: dict[str, int], linear: bool, quantity_source: str) -> dict:
    """Provenance recorded in ``AnnData.uns`` for one modality."""
    return {
        "missing_value_convention": LINEAR_CONVENTION if linear else LOG2_CONVENTION,
        "n_source_zero_to_nan": stats["n_source_zero_to_nan"],
        "n_unstored_to_nan": stats["n_unstored_to_nan"] + stats.get("n_source_nan", 0),
        "quantity_source": quantity_source,
    }


# ---------------------------------------------------------------------------
# Guards G1-G11 (G12-G14 are TMT-specific, in tmt_to_anndata.py)
# ---------------------------------------------------------------------------


def check_failed_modalities(uns) -> None:
    """G1: qpx >= 1.1.3 records modality build errors; the record must exist and be empty."""
    if "qpx_failed_modalities" not in uns:
        raise GuardError("G1", "uns['qpx_failed_modalities'] is missing (qpx older than 1.1.3?); cannot verify the build")
    failed = dict(uns["qpx_failed_modalities"])
    if failed:
        raise GuardError("G1", f"MuData modalities failed to build: {failed}")


def check_modality_present(modalities, name: str, guard: str, hint: str) -> None:
    """G2 (proteins) / G7 (precursors): a required modality exists."""
    if name not in modalities:
        raise GuardError(guard, f"no '{name}' modality (built: {sorted(modalities)}); {hint}")


def check_nonempty(name: str, dense: np.ndarray) -> None:
    """G4: at least one obs, one var and one non-NaN value."""
    if dense.shape[0] == 0 or dense.shape[1] == 0 or not np.any(~np.isnan(dense)):
        raise GuardError("G4", f"'{name}' matrix is empty: shape {dense.shape}, {int((~np.isnan(dense)).sum())} values")


def check_sample_accession(name: str, obs: pd.DataFrame) -> None:
    """G5: every obs row carries an SDRF ``sample_accession``."""
    if "sample_accession" not in obs.columns:
        raise GuardError("G5", f"'{name}' obs has no sample_accession column: no SDRF sample matched any run")
    empty = sorted(obs.index[is_empty(obs["sample_accession"])].astype(str))
    if empty:
        raise GuardError("G5", f"'{name}' obs has empty sample_accession for {empty}: runs not in the SDRF")


def is_empty(values: pd.Series) -> pd.Series:
    """True where a metadata value is missing: None, pd.NA, NaN, or ''/'None'/'<NA>'/'nan' after strip."""
    as_text = values.astype(str).str.strip().str.lower()
    return values.isna() | as_text.isin(EMPTY_TOKENS)


def check_obs_names(name: str, obs_names, expected) -> None:
    """G6: obs_names equal the expected run/experiment set."""
    got, want = set(map(str, obs_names)), set(map(str, expected))
    if got != want:
        raise GuardError(
            "G6",
            f"'{name}' obs_names differ from the expected set: missing={sorted(want - got)}, "
            f"unexpected={sorted(got - want)}",
        )


def check_nan_count(name: str, dense: np.ndarray, expected_nan: int) -> None:
    """G8: the NaN count equals the missing count derived from the source."""
    n_nan = int(np.isnan(dense).sum())
    if n_nan != expected_nan:
        raise GuardError("G8", f"'{name}' has {n_nan} NaN but the source has {expected_nan} missing cells")


def check_value_range(name: str, dense: np.ndarray, linear: bool) -> None:
    """G9: no inf anywhere; no value <= 0 on a linear scale.

    Run on the FINAL matrices (X and every layer, after collapse / reorder /
    pivot), right before writing: ``densify_with_nan`` already maps <= 0 to NaN,
    so G9 only guards against a later step re-introducing 0 or inf.
    """
    n_inf = int(np.isinf(dense).sum())
    if n_inf:
        raise GuardError("G9", f"'{name}' contains {n_inf} infinite values")
    if linear:
        n_nonpos = int((dense <= 0).sum())
        if n_nonpos:
            raise GuardError("G9", f"'{name}' contains {n_nonpos} values <= 0 on a linear scale (zero leak)")


def check_aligned_obs(obs_by_modality: dict[str, pd.Index]) -> None:
    """G10: every modality has the same obs_names in the same order."""
    items = list(obs_by_modality.items())
    ref_name, ref = items[0]
    for name, names in items[1:]:
        if list(map(str, names)) != list(map(str, ref)):
            raise GuardError(
                "G10",
                f"'{name}' obs_names {list(names)} are not aligned with '{ref_name}' obs_names {list(ref)}",
            )


def check_fraction_rows_identical(dense: np.ndarray, obs_names, exp_to_runs: dict[str, list[str]]) -> None:
    """G11: all fraction rows of one experiment are identical (NaN-aware) before the collapse."""
    position = {str(n): i for i, n in enumerate(obs_names)}
    for experiment, runs in exp_to_runs.items():
        rows = dense[[position[r] for r in runs]]
        if not all(np.array_equal(rows[0], row, equal_nan=True) for row in rows[1:]):
            raise GuardError(
                "G11",
                f"fraction rows of experiment '{experiment}' ({runs}) differ; refusing to collapse them to one obs",
            )


# ---------------------------------------------------------------------------
# FragPipe helpers
# ---------------------------------------------------------------------------


def read_tsv_header(path) -> list[str]:
    """Return the header row of a TSV file."""
    with open(path, newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def fragpipe_header_experiments(columns) -> list[str]:
    """Experiments named by ``<exp> Spectral Count`` columns of combined_protein.tsv."""
    suffix = " Spectral Count"
    return sorted(
        c[: -len(suffix)]
        for c in columns
        if c.endswith(suffix) and not c.endswith((" Unique Spectral Count", " Total Spectral Count"))
        and not c.startswith("Combined ")
    )


def read_experiment_runs(sdrf_path, annotation_path=None) -> dict[str, list[str]]:
    """Map each FragPipe experiment to its run stems.

    With FragPipe's ``experiment_annotation.tsv`` (``file``/``sample`` columns, the
    IonQuant experiment names ``<experiment>_<bioreplicate>``) that file is
    authoritative, as in qpx. Otherwise the SDRF ``source name`` is the
    experiment (pipeline naming) and ``comment[data file]`` gives its runs; rows
    sharing a source name are runs of one experiment. Run order follows the file.
    """
    if annotation_path:
        table = pd.read_csv(annotation_path, sep="\t", dtype=str, keep_default_na=False)
        exp_col, file_col = "sample", "file"
        source = f"experiment annotation {Path(annotation_path).name}"
    else:
        table = pd.read_csv(sdrf_path, sep="\t", dtype=str, keep_default_na=False)
        table.columns = [c.strip().lower() for c in table.columns]
        exp_col, file_col = "source name", "comment[data file]"
        source = f"SDRF {Path(sdrf_path).name}"
    missing = {exp_col, file_col} - set(table.columns)
    if missing:
        raise GuardError("G3", f"{source} lacks column(s) {sorted(missing)}")
    mapping: dict[str, list[str]] = {}
    for experiment, data_file in zip(table[exp_col].str.strip(), table[file_col].str.strip()):
        if not experiment or not data_file:
            raise GuardError("G3", f"{source} has an empty '{exp_col}' or '{file_col}' value")
        runs = mapping.setdefault(experiment, [])
        stem = run_file_stem(data_file)
        if stem not in runs:
            runs.append(stem)
    return mapping


def sdrf_run_stems(sdrf_path) -> list[str]:
    """Run stems of every SDRF ``comment[data file]`` (the expected run set)."""
    table = pd.read_csv(sdrf_path, sep="\t", dtype=str, keep_default_na=False)
    table.columns = [c.strip().lower() for c in table.columns]
    if "comment[data file]" not in table.columns:
        raise GuardError("G6", f"SDRF {Path(sdrf_path).name} has no comment[data file] column")
    return list(dict.fromkeys(run_file_stem(v) for v in table["comment[data file]"]))


def precheck_fragpipe_columns(columns, exp_to_runs: dict[str, list[str]]) -> str:
    """G3: every expected experiment has a MaxLFQ or Total intensity column.

    Runs before ``qpxc``, which silently drops the whole protein view when these
    columns are absent. Returns the quantity-source description for ``uns``.
    """
    header_experiments = fragpipe_header_experiments(columns)
    expected = sorted(exp_to_runs)
    unmatched = sorted(set(expected) - set(header_experiments))
    if unmatched:
        hint = ""
        if any(f"{e}_1" in header_experiments for e in unmatched):
            hint = "; IonQuant naming '<experiment>_<bioreplicate>' needs --experiment-annotation-file"
        raise GuardError(
            "G3",
            f"experiments {unmatched} have no columns in combined_protein.tsv "
            f"(experiments found in the header: {header_experiments}){hint}",
        )
    extra = sorted(set(header_experiments) - set(expected))
    if extra:
        raise GuardError("G3", f"combined_protein.tsv experiments {extra} are not in the SDRF/annotation")
    column_set = set(columns)
    lacking = [
        e for e in expected if f"{e} MaxLFQ Intensity" not in column_set and f"{e} Total Intensity" not in column_set
    ]
    if lacking:
        raise GuardError(
            "G3",
            f"combined_protein.tsv has no MaxLFQ Intensity / Total Intensity column for experiments {lacking} "
            f"(only '<exp> Intensity'); qpx would silently drop the proteins modality. {FRAGPIPE_MAXLFQ_HINT}",
        )
    if all(f"{e} Total Intensity" in column_set for e in expected):
        return "<exp> Total Intensity (combined_protein.tsv; qpx uses <exp> MaxLFQ Intensity where Total <= 0)"
    return "<exp> MaxLFQ Intensity (combined_protein.tsv)"


def collapse_fragpipe_obs(dense: np.ndarray, obs: pd.DataFrame, exp_to_runs: dict[str, list[str]]):
    """Collapse fraction rows to one row per experiment. Call only after G11 passed.

    obs columns that agree across the fractions keep their value; columns that
    differ (run_accession, fraction, ...) are joined with ``;``. ``run_file_names``
    lists the member runs. Returns ``(dense, obs)`` indexed by experiment.
    """
    position = {str(n): i for i, n in enumerate(obs.index)}
    experiments = sorted(exp_to_runs)
    rows, records = [], []
    for experiment in experiments:
        runs = exp_to_runs[experiment]
        rows.append(dense[position[runs[0]]])
        member = obs.iloc[[position[r] for r in runs]]
        record = {col: _join_unique(member[col]) for col in obs.columns}
        record["run_file_names"] = ";".join(runs)
        records.append(record)
    collapsed_obs = pd.DataFrame(records, index=pd.Index(experiments, name="experiment"))
    return np.vstack(rows).astype(np.float32), collapsed_obs


def _join_unique(values: pd.Series):
    unique = list(dict.fromkeys(values.tolist()))
    return unique[0] if len(unique) == 1 else ";".join(map(str, unique))


def fragpipe_protein_source_missing(pg_path, experiments, var_names) -> dict[str, int]:
    """Count the missing protein cells in combined_protein.tsv itself (FragPipe G8 source).

    qpx drops FragPipe cells <= 0 before the matrix is built, so the matrix alone
    cannot tell a source 0 from an absent cell. Over the rows kept as ``var``
    (anchor = first ``Protein`` token, as in qpx) and the qpx primary columns
    (``<exp> Total Intensity``, else ``<exp> MaxLFQ Intensity``), a cell is a
    source zero when no column is > 0 and one holds a number <= 0, and absent
    when every column is NA/empty. One cell per experiment and protein.
    """
    table = pd.read_csv(pg_path, sep="\t", dtype=str, keep_default_na=False)
    anchors = table["Protein"].str.split(",").str[0].str.strip()
    rows = table[anchors.isin(set(map(str, var_names)))]
    if len(rows) != len(var_names):
        raise GuardError(
            "G8", f"combined_protein.tsv has {len(rows)} rows for the {len(var_names)} protein vars (duplicate anchors?)"
        )
    n_zero = n_absent = 0
    for experiment in experiments:
        columns = [c for c in (f"{experiment} Total Intensity", f"{experiment} MaxLFQ Intensity") if c in rows.columns]
        values = rows[columns].apply(lambda col: pd.to_numeric(col.str.strip().replace({"": None, "NA": None})))
        positive = (values > 0).any(axis=1)
        absent = values.isna().all(axis=1)
        n_zero += int((~positive & ~absent).sum())
        n_absent += int(absent.sum())
    return {"n_source_zero_to_nan": n_zero, "n_unstored_to_nan": n_absent}


def align_fragpipe_precursor_obs(precursor_obs_names, protein_obs: pd.DataFrame, exp_to_runs, collapsed: bool):
    """Return ``(order, obs)`` that aligns precursors (keyed by experiment) with proteins.

    Not collapsed: each experiment has one run, so precursor rows are relabelled
    to that run stem. Collapsed: proteins are already keyed by experiment.
    ``order`` holds the precursor row positions in protein obs order; ``obs`` is
    the protein obs (SDRF metadata) for the precursors modality.
    """
    if collapsed:
        label = {e: e for e in exp_to_runs}
    else:
        label = {e: runs[0] for e, runs in exp_to_runs.items()}
    relabelled = [label[str(e)] for e in precursor_obs_names]
    position = {name: i for i, name in enumerate(relabelled)}
    order = [position[str(name)] for name in protein_obs.index]
    return order, protein_obs.copy()


def fragpipe_dropped_proteins(pg_path, experiments, var_names) -> dict[str, int]:
    """Account for combined_protein.tsv rows that qpx did not keep as protein vars (D4).

    qpx drops a protein whose MaxLFQ (or Total) intensity is <= 0 / empty in every
    experiment. Fails (G8) if a dropped row has any positive value, which would be
    real data lost. Also counts the cells that are NaN in X although the
    ``<exp> Intensity`` column is > 0 (MaxLFQ could not be estimated there).
    """
    table = pd.read_csv(pg_path, sep="\t", dtype=str, keep_default_na=False)
    anchors = table["Protein"].str.split(",").str[0].str.strip()
    kept = anchors.isin(set(map(str, var_names)))

    def numeric(columns):
        values = table[[c for c in columns if c in table.columns]]
        return values.apply(lambda col: pd.to_numeric(col.str.strip().replace({"": None, "NA": None})))

    primary = numeric([c for e in experiments for c in (f"{e} Total Intensity", f"{e} MaxLFQ Intensity")])
    lost = sorted(anchors[~kept & (primary > 0).any(axis=1)])
    if lost:
        raise GuardError("G8", f"qpx dropped proteins with a positive MaxLFQ/Total intensity: {lost[:5]}")
    n_nan_with_intensity = 0
    for experiment in experiments:
        quantified = numeric([f"{experiment} Total Intensity", f"{experiment} MaxLFQ Intensity"])
        intensity = numeric([f"{experiment} Intensity"])
        missing = ~(quantified > 0).any(axis=1)
        n_nan_with_intensity += int((kept & missing & (intensity > 0).any(axis=1)).sum())
    raw_intensity = numeric([f"{e} Intensity" for e in experiments])
    return {
        "n_proteins_all_missing_dropped": int((~kept).sum()),
        "n_dropped_with_intensity": int((~kept & (raw_intensity > 0).any(axis=1)).sum()),
        "n_nan_with_intensity": n_nan_with_intensity,
    }


# ---------------------------------------------------------------------------
# obs metadata (SDRF / annotation join)
# ---------------------------------------------------------------------------


def sanitize_column(name: str) -> str:
    """Deterministic obs column name for an SDRF header.

    ``characteristics[organism part]`` -> ``organism_part``, ``factor value[x]`` ->
    ``factor_x``, ``comment[x]`` -> ``comment_x``, anything else lower-cased;
    non-alphanumerics become ``_``. pandas' ``.1`` suffix on a repeated header
    becomes ``_1``.
    """
    text = name.strip().lower()
    for prefix, replacement in (("characteristics[", ""), ("factor value[", "factor_"), ("comment[", "comment_")):
        if text.startswith(prefix):
            text = replacement + text[len(prefix):].replace("]", "", 1)
            break
    return re.sub(r"[^0-9a-z]+", "_", text).strip("_")


def sdrf_obs_table(sdrf_path) -> pd.DataFrame:
    """Every SDRF column, sanitized, one row per run stem (``comment[data file]``)."""
    table = pd.read_csv(sdrf_path, sep="\t", dtype=str, keep_default_na=False)
    data_file = [c for c in table.columns if c.strip().lower() == "comment[data file]"]
    if len(data_file) != 1:
        raise GuardError("G5", f"SDRF {Path(sdrf_path).name} needs exactly one comment[data file] column")
    stems = table[data_file[0]].map(run_file_stem)
    repeated = sorted(stems[stems.duplicated()].unique())
    if repeated:
        raise GuardError("G5", f"SDRF {Path(sdrf_path).name} lists runs {repeated} on more than one row")
    names = [sanitize_column(c) for c in table.columns]
    clashes = sorted({n for n in names if names.count(n) > 1})
    if clashes:
        raise GuardError("G5", f"SDRF {Path(sdrf_path).name} columns collide after sanitizing: {clashes}")
    table.columns = names
    return table.apply(lambda col: col.str.strip()).set_index(stems.rename(None))


def annotation_obs_table(annotation_path) -> pd.DataFrame:
    """FragPipe experiment_annotation.tsv ``condition`` / ``replicate`` per run stem."""
    table = pd.read_csv(annotation_path, sep="\t", dtype=str, keep_default_na=False)
    columns = [c for c in ("condition", "replicate") if c in table.columns]
    stems = table["file"].map(run_file_stem)
    return table[columns].apply(lambda col: col.str.strip()).set_index(stems.rename(None))


def join_obs_metadata(obs: pd.DataFrame, extra: pd.DataFrame, clash_prefix: str) -> pd.DataFrame:
    """Left-join *extra* (indexed by run stem) onto *obs* (indexed by run stem).

    A column whose name is already in *obs* is added as ``<clash_prefix><name>``.
    Every obs row must have a metadata row (G5); values are strings ('' if empty).
    """
    unmatched = sorted(set(map(str, obs.index)) - set(map(str, extra.index)))
    if unmatched:
        raise GuardError("G5", f"obs {unmatched} have no metadata row (runs not in the SDRF/annotation)")
    renamed = extra.rename(columns={c: f"{clash_prefix}{c}" for c in extra.columns if c in obs.columns})
    joined = obs.join(renamed.reindex(obs.index.astype(str)).set_axis(obs.index))
    return joined.astype({c: str for c in renamed.columns})


# ---------------------------------------------------------------------------
# Protein identifiers (var columns shared by every engine)
# ---------------------------------------------------------------------------


def protein_accession(token: str) -> str:
    """Bare accession: ``sp|P12345|NAME`` -> ``P12345``; anything else unchanged (stripped)."""
    parts = str(token).strip().split("|")
    return parts[1] if len(parts) >= 3 else str(token).strip()


def is_contaminant(token: str) -> bool:
    """True for a contaminant entry, by :data:`CONTAMINANT_PREFIXES`."""
    return str(token).strip().lower().startswith(CONTAMINANT_PREFIXES)


def protein_var_columns(leaders, members) -> pd.DataFrame:
    """``protein_accession``, ``protein_group`` (``;``-joined accessions) and ``is_contaminant``.

    *leaders* is the leading protein per var; *members* the list of every protein
    of the group (leader first). ``is_contaminant`` is true when any member is one.
    """
    return pd.DataFrame(
        {
            "protein_accession": [protein_accession(leader) for leader in leaders],
            "protein_group": [";".join(dict.fromkeys(protein_accession(m) for m in group)) for group in members],
            "is_contaminant": [any(is_contaminant(m) for m in group) for group in members],
        }
    )


def split_members(value: str, separator: str) -> list[str]:
    return [part.strip() for part in str(value).split(separator) if part.strip()]


# ---------------------------------------------------------------------------
# DIA-NN pg_matrix
# ---------------------------------------------------------------------------

PG_MATRIX_META = ("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description")


def canonical_membership(accessions) -> str:
    """Order-independent protein-group key, as qpx ``_canonical_pg_membership``."""
    return ";".join(sorted({a.strip() for a in accessions if a and a.strip()}))


def read_pg_matrix(path) -> tuple[pd.DataFrame, pd.Series]:
    """DIA-NN pg_matrix values indexed by canonical membership, columns = run stems.

    Returns ``(values, protein_group)``; ``values`` is float64 with NaN for empty
    cells, ``protein_group`` is DIA-NN's own ``Protein.Group`` string.
    """
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "Protein.Group" not in table.columns:
        raise GuardError("G17", f"{Path(path).name} has no Protein.Group column; not a DIA-NN pg_matrix")
    keys = table["Protein.Group"].map(lambda v: canonical_membership(v.split(";")))
    repeated = sorted(keys[keys.duplicated()].unique())
    if repeated:
        raise GuardError("G17", f"{Path(path).name} repeats protein groups {repeated[:5]}")
    run_columns = [c for c in table.columns if c not in PG_MATRIX_META]
    values = table[run_columns].apply(lambda col: pd.to_numeric(col.str.strip().replace({"": None})))
    values.columns = [run_file_stem(c) for c in run_columns]
    values.index = keys.rename(None)
    return values, table["Protein.Group"].set_axis(keys.rename(None))


def check_matches_pg_matrix(dense: np.ndarray, obs_names, memberships, pg_values: pd.DataFrame) -> None:
    """G17: proteins X equals pg_matrix cell for cell (same non-NaN cells, float32 values)."""
    expected = pg_values.reindex(index=list(memberships), columns=list(map(str, obs_names))).to_numpy().T
    expected = np.where(expected > 0, expected, np.nan).astype(np.float32)
    n_x, n_pg = int((~np.isnan(dense)).sum()), int((~np.isnan(expected)).sum())
    n_pg_matrix = int((pg_values.to_numpy() > 0).sum())
    if n_x != n_pg or n_pg != n_pg_matrix:
        raise GuardError(
            "G17", f"proteins X has {n_x} values but pg_matrix has {n_pg_matrix} ({n_pg} on the exported axes)"
        )
    if not np.allclose(dense, expected, rtol=1e-6, atol=0, equal_nan=True):
        raise GuardError("G17", "proteins X values differ from pg_matrix beyond float32 tolerance")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------


def summary_rows(prefix: str, name: str, dense: np.ndarray, uns: dict, layers=None) -> dict[str, str]:
    """Summary rows for one modality: ``count.*`` rows are data-dependent, the rest are stable.

    Every integer ``uns`` entry named ``n_*`` becomes a ``count.`` row; every
    ``uns`` string entry in :data:`SUMMARY_TEXT_KEYS` a stable row. Each layer in
    *layers* (``{name: matrix}``) adds its source (``uns['layers'][name]``) and counts.
    """
    key = f"{prefix}{name}"
    rows = {f"{key}.{k}": str(uns[k]) for k in SUMMARY_TEXT_KEYS if k in uns}
    rows.update({f"count.{key}.{k}": str(v) for k, v in uns.items() if k.startswith("n_") and isinstance(v, numbers.Integral)})
    rows.update(_matrix_counts(f"count.{key}", dense))
    for layer, matrix in (layers or {}).items():
        rows[f"{key}.layers.{layer}.quantity_source"] = uns["layers"][layer]
        rows.update(_matrix_counts(f"count.{key}.layers.{layer}", matrix))
    return rows


SUMMARY_TEXT_KEYS = ("missing_value_convention", "quantity_source", "scale", "scale_rule")


def _matrix_counts(key: str, dense: np.ndarray) -> dict[str, str]:
    return {
        f"{key}.n_obs": str(dense.shape[0]),
        f"{key}.n_vars": str(dense.shape[1]),
        f"{key}.n_nan": str(int(np.isnan(dense).sum())),
        f"{key}.n_non_nan": str(int((~np.isnan(dense)).sum())),
    }


def write_summary_tsv(path, rows: dict[str, str]) -> None:
    """Write ``key<TAB>value`` rows sorted by key (deterministic)."""
    for key, value in rows.items():
        if "\t" in key or "\n" in key or "\t" in str(value) or "\n" in str(value):
            raise ValueError(f"summary row {key!r} contains a tab or newline")
    with open(path, "w") as handle:
        for key in sorted(rows):
            handle.write(f"{key}\t{rows[key]}\n")
