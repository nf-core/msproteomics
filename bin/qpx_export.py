#!/usr/bin/env python3
"""Build guarded .h5mu / .h5ad files from a qpx (qpxc convert) output folder.

Runs after ``qpxc convert diann`` or ``qpxc convert fragpipe``. It assembles the
MuData with ``qpx.mudata.build_mudata``, runs the silent-failure guards
(see ``qpx_export_lib.py``), converts every matrix to dense float32 with NaN
for missing values, and only then writes:

    <prefix>.h5mu            all modalities (precursors, proteins)
    <prefix>.proteins.h5ad   the proteins modality
    <prefix>.summary.tsv     sorted key<TAB>value rows (count.* rows are data-dependent)

DIA-NN values (user decision 2026-09-27, docs/plans/2026-09-27_h5ad_decisions.md "DIA X"):
``proteins.X`` = PG.MaxLFQ read from DIA-NN's own ``pg_matrix.tsv`` (``--pg-matrix-file``;
DIA-NN applies its precursor ``Global.Q.Value`` and protein-group FDR filter there),
``layers['raw']`` = PG.Quantity from qpx on the same cells; ``precursors.X`` =
Precursor.Normalised on DIA-NN's pr_matrix cells (Global.Q.Value and Global.PG.Q.Value
<= ``--matrix-qvalue``; G18 against ``--pr-matrix-file`` when given), ``layers['raw']`` =
Precursor.Quantity on the same cells.
obs carries every SDRF column (sanitized names) and, for allinone FragPipe, the
annotation ``condition`` / ``replicate``.

``--precheck`` (FragPipe) only runs guard G3 on combined_protein.tsv, before
``qpxc`` is started, because qpxc exits 0 while dropping the protein view when
the MaxLFQ / Total intensity columns are missing.

Any guard failure prints ``ERROR: G<n>: ...`` and exits 3 without writing.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from qpx_export_lib import (
    GUARD_EXIT_CODE,
    GuardError,
    align_fragpipe_precursor_obs,
    annotation_obs_table,
    canonical_membership,
    check_aligned_obs,
    check_failed_modalities,
    check_fraction_rows_identical,
    check_matches_pg_matrix,
    check_modality_present,
    check_nan_count,
    check_nonempty,
    check_obs_names,
    check_sample_accession,
    check_value_range,
    collapse_fragpipe_obs,
    densify_with_nan,
    fragpipe_dropped_proteins,
    fragpipe_protein_source_missing,
    join_obs_metadata,
    modality_uns,
    precheck_fragpipe_columns,
    protein_var_columns,
    read_experiment_runs,
    read_pg_matrix,
    read_tsv_header,
    run_file_stem,
    sdrf_obs_table,
    sdrf_run_stems,
    split_members,
    summary_rows,
    write_summary_tsv,
)

MODALITIES = ("precursors", "proteins")
GUARDS_RUN = {
    "diann": "G1,G2,G4,G5,G6,G7,G8,G9,G10,G16,G17",
    "fragpipe": "G1,G2,G3,G4,G5,G6,G7,G8,G9,G10,G11,G16",
}

# Source rows, distinct (run, var) cells, and distinct cells with a positive value
# in the qpx Parquet views, keyed like build_mudata pivots them. rows != cells means
# qpx would sum duplicate cells (G16); positive cells are the G8 source count.
CELL_COUNT_QUERIES = {
    "proteins": """
        SELECT COUNT(*), COUNT(DISTINCT k), COUNT(DISTINCT CASE WHEN v > 0 THEN k END) FROM (
            SELECT gr || chr(31) || pg.anchor_protein AS k, pg.intensity AS v
            FROM read_parquet($path) AS pg, UNNEST(list_distinct(pg.grouped_runs)) AS _g(gr)
            WHERE pg.label IS NOT NULL)
    """,
    "precursors": """
        SELECT COUNT(*), COUNT(DISTINCT k), COUNT(DISTINCT CASE WHEN v > 0 THEN k END) FROM (
            SELECT f.run_file_name || chr(31) || f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS k,
                   i.intensity AS v
            FROM read_parquet($path) AS f, UNNEST(f.intensities) AS _t(i))
    """,
}
PARQUET_VIEW = {"proteins": "pg", "precursors": "feature"}
DIANN_PG_METHOD_QUERY = """
    SELECT DISTINCT c.cv_value
    FROM read_parquet($path) AS pg, UNNEST(pg.cv_params) AS _t(c)
    WHERE c.cv_name = 'quantification_method'
"""
DIANN_MEMBERSHIP_QUERY = "SELECT DISTINCT anchor_protein, pg_accessions FROM read_parquet($path)"
# DIA-NN pr_matrix = report rows with Global.Q.Value <= q AND Global.PG.Q.Value <= q ("1% precursor and
# protein group FDR", diann_report.log.txt), values = Precursor.Normalised. Reproduced cell for cell on
# EXP26000884 (104,050 cells / 26,948 precursors). qpx keeps Global.Q.Value as the additional score
# 'global_qvalue' and Global.PG.Q.Value as pg_global_qvalue.
DIANN_NORMALISED_QUERY = """
    SELECT f.run_file_name, f.peptidoform || '|' || CAST(f.charge AS VARCHAR) AS precursor_id,
           list_filter(a.intensities, x -> x.intensity_name = 'precursor_normalised')[1].intensity_value AS value,
           list_filter(f.additional_scores, x -> x.score_name = 'global_qvalue')[1].score_value AS global_qvalue,
           f.pg_global_qvalue
    FROM read_parquet($path) AS f, UNNEST(f.additional_intensities) AS _a(a)
"""
DIANN_SOURCES = {
    "proteins": "PG.MaxLFQ (DIA-NN pg_matrix; normalised)",
    "proteins.raw": "PG.Quantity (DIA-NN report via qpx; not normalised; pg_matrix cells only)",
    "precursors": "Precursor.Normalised (DIA-NN report via qpx, pr_matrix filter "
    "Global.Q.Value & Global.PG.Q.Value <= {q}; normalised)",
    "precursors.raw": "Precursor.Quantity (DIA-NN report via qpx; not normalised; pr_matrix cells only)",
}
MISSING_MODALITY_HINT = {
    "diann": {
        "proteins": "DIA-NN pg_matrix missing (--pg-matrix-path) or empty after the q-value filter",
        "precursors": "DIA-NN report produced no precursor features",
    },
    "fragpipe": {
        "proteins": "combined_protein.tsv has no MaxLFQ/Total intensity columns or no protein passed",
        "precursors": "combined_ion.tsv (--ion-file) produced no features",
    },
}


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--engine", required=True, choices=["diann", "fragpipe"], help="qpxc converter that ran")
    parser.add_argument("--sdrf-file", required=True, help="SDRF passed to qpxc (expected run set, obs metadata)")
    parser.add_argument("--qpx-folder", help="qpxc --output-folder")
    parser.add_argument("--prefix", help="qpxc --output-prefix; also the output file prefix")
    parser.add_argument("--pg-matrix-file", help="DIA-NN pg_matrix.tsv (required for --engine diann)")
    parser.add_argument("--pr-matrix-file", help="DIA-NN pr_matrix.tsv (optional): guard G18 compares precursors.X to it")
    parser.add_argument("--matrix-qvalue", type=float, default=0.01,
                        help="DIA-NN matrix FDR (--qvalue of the DIA-NN run; default 0.01) for the precursor filter")
    parser.add_argument("--pg-file", help="FragPipe combined_protein.tsv (required for --engine fragpipe)")
    parser.add_argument("--experiment-annotation-file", help="FragPipe experiment_annotation.tsv (allinone)")
    parser.add_argument(
        "--sample-map",
        help="FragPipe pipeline mode: headerless '<file><TAB><group>' TSV (buildFileExperimentMap); "
        "authoritative run -> experiment mapping",
    )
    parser.add_argument("--outdir", default=".", help="output directory (default: .)")
    parser.add_argument("--precheck", action="store_true", help="FragPipe only: run guard G3 and exit")
    args = parser.parse_args(argv)
    if args.engine == "fragpipe" and not args.pg_file:
        parser.error("--engine fragpipe requires --pg-file")
    if args.sample_map and args.experiment_annotation_file:
        parser.error("--sample-map and --experiment-annotation-file are mutually exclusive")
    if args.sample_map and args.engine != "fragpipe":
        parser.error("--sample-map is only defined for --engine fragpipe")
    if args.precheck and args.engine != "fragpipe":
        parser.error("--precheck is only defined for --engine fragpipe")
    if not args.precheck and not (args.qpx_folder and args.prefix):
        parser.error("--qpx-folder and --prefix are required unless --precheck is given")
    if args.engine == "diann" and not args.pg_matrix_file:
        parser.error("--engine diann requires --pg-matrix-file")
    return args


def experiment_runs(args) -> dict[str, list[str]]:
    return read_experiment_runs(args.sdrf_file, args.experiment_annotation_file, args.sample_map)


def run_precheck(args) -> str:
    """G3 on combined_protein.tsv against the sample-map / annotation / SDRF experiments."""
    return precheck_fragpipe_columns(read_tsv_header(args.pg_file), experiment_runs(args))


def query_rows(sql: str, path: Path) -> list[tuple]:
    import duckdb

    with duckdb.connect() as con:
        return con.execute(sql, {"path": str(path)}).fetchall()


def parquet_path(args, modality: str) -> Path:
    path = Path(args.qpx_folder) / f"{args.prefix}.{PARQUET_VIEW[modality]}.parquet"
    if not path.is_file():
        raise GuardError("G8", f"source view {path} is missing; cannot count source cells for '{modality}'")
    return path


def build_qpx_mudata(folder: str):
    from qpx.dataset import Dataset
    from qpx.mudata import build_mudata

    dataset = Dataset(folder)
    mdata = build_mudata(dataset)
    dataset.close()
    return mdata


def densify_modalities(args, mdata) -> dict[str, tuple[np.ndarray, dict]]:
    """G16 + G8 per modality: no duplicate source cells; NaN count matches the Parquet source."""
    result = {}
    for name in MODALITIES:
        n_rows, n_cells, n_positive = query_rows(CELL_COUNT_QUERIES[name], parquet_path(args, name))[0]
        if n_rows != n_cells:
            raise GuardError(
                "G16",
                f"{PARQUET_VIEW[name]}.parquet has {n_rows} rows for {n_cells} distinct (run, {name[:-1]}) cells; "
                "qpx would sum the duplicates",
            )
        dense, stats = densify_with_nan(mdata.mod[name].X, linear=True)
        check_nan_count(name, dense, dense.size - int(n_positive))
        result[name] = (dense, stats)
    return result


def pivot_long(frame: pd.DataFrame, obs_names, var_names, name: str) -> np.ndarray:
    """(obs, var, value) rows -> dense float64 matrix on the given axes (NaN = no row); G16 on duplicates."""
    obs_col, var_col, value_col = frame.columns
    duplicated = frame.duplicated(subset=[obs_col, var_col])
    if duplicated.any():
        raise GuardError("G16", f"'{name}' has {int(duplicated.sum())} duplicate (run, var) source rows")
    rows = pd.Index(obs_names).get_indexer(frame[obs_col])
    cols = pd.Index(var_names).get_indexer(frame[var_col])
    if (rows < 0).any() or (cols < 0).any():
        raise GuardError("G16", f"'{name}' has {int(((rows < 0) | (cols < 0)).sum())} source rows off the qpx axes")
    dense = np.full((len(obs_names), len(var_names)), np.nan)
    dense[rows, cols] = frame[value_col].to_numpy(dtype=np.float64)
    return dense


def diann_proteins(args, proteins, raw: np.ndarray):
    """proteins.X from pg_matrix (G17), layers['raw'] = qpx PG.Quantity on the pg_matrix cells."""
    methods = sorted(row[0] for row in query_rows(DIANN_PG_METHOD_QUERY, parquet_path(args, "proteins")))
    if methods != ["PG.Quantity"]:
        raise GuardError("G2", f"pg.parquet quantification_method is {methods}, expected ['PG.Quantity'] for layers['raw']")
    pg_values, pg_group = read_pg_matrix(args.pg_matrix_file)
    check_obs_names("pg_matrix runs", pg_values.columns, proteins.obs_names)
    key_of: dict[str, str] = {}
    for anchor, accessions in query_rows(DIANN_MEMBERSHIP_QUERY, parquet_path(args, "proteins")):
        key = canonical_membership(accessions)
        if key_of.setdefault(anchor, key) != key:
            raise GuardError("G16", f"qpx anchor '{anchor}' stands for two protein groups ({key_of[anchor]}, {key})")
    keys = [key_of[str(a)] for a in proteins.var_names]
    if len(set(keys)) != len(keys):
        raise GuardError("G16", "two qpx protein vars share one protein-group membership")
    unmapped = sorted(set(pg_values.index) - set(keys))
    if unmapped:
        raise GuardError("G17", f"pg_matrix protein groups {unmapped[:5]} ({len(unmapped)}) are not in the qpx proteins view")
    keep = [i for i, key in enumerate(keys) if key in pg_values.index]
    kept_keys = [keys[i] for i in keep]
    values = pg_values.reindex(index=kept_keys, columns=list(map(str, proteins.obs_names))).to_numpy().T
    dense, stats = densify_with_nan(values, linear=True)
    check_matches_pg_matrix(dense, proteins.obs_names, kept_keys, pg_values)
    raw = raw[:, keep]
    lacking = int((~np.isnan(dense) & np.isnan(raw)).sum())
    if lacking:
        raise GuardError("G17", f"{lacking} pg_matrix cells have no PG.Quantity in the qpx proteins view")
    n_masked = int((np.isnan(dense) & ~np.isnan(raw)).sum())
    raw[np.isnan(dense)] = np.nan
    var = proteins.var.iloc[keep].copy()
    groups = [split_members(pg_group[key], ";") for key in kept_keys]
    ids = protein_var_columns(list(var.index), groups).set_index(var.index)
    var = var.join(ids)
    uns = modality_uns(stats, linear=True, quantity_source=DIANN_SOURCES["proteins"])
    uns.update({
        "layers": {"raw": DIANN_SOURCES["proteins.raw"]},
        "n_report_only_groups_dropped": len(keys) - len(keep),
        "n_raw_cells_outside_pg_matrix": n_masked,
    })
    return dense, {"raw": raw}, var, uns


def diann_precursors(args, precursors, raw: np.ndarray):
    """precursors.X = Precursor.Normalised on DIA-NN's pr_matrix cells; layers['raw'] = Precursor.Quantity there."""
    frame = pd.DataFrame(
        query_rows(DIANN_NORMALISED_QUERY, parquet_path(args, "precursors")),
        columns=["run", "precursor", "value", "global_qvalue", "pg_global_qvalue"],
    )
    for column in ("value", "global_qvalue", "pg_global_qvalue"):
        if frame[column].isna().any():
            raise GuardError("G7", f"feature.parquet lacks '{column}' for {int(frame[column].isna().sum())} precursor rows")
    passed = frame[(frame["global_qvalue"] <= args.matrix_qvalue) & (frame["pg_global_qvalue"] <= args.matrix_qvalue)]
    if passed.empty:
        raise GuardError("G7", f"no precursor passes Global.Q.Value and Global.PG.Q.Value <= {args.matrix_qvalue}")
    keep = pd.Index(precursors.var_names).isin(passed["precursor"])
    var_names = precursors.var_names[keep]
    values = pivot_long(passed[["run", "precursor", "value"]], precursors.obs_names, var_names, "precursors")
    dense, stats = densify_with_nan(values, linear=True)
    check_nan_count("precursors", dense, dense.size - int((passed["value"] > 0).sum()))
    identified = ~np.isnan(pivot_long(passed[["run", "precursor", "global_qvalue"]], precursors.obs_names, var_names, "precursors"))
    raw = raw[:, keep]
    n_masked = int((~identified & ~np.isnan(raw)).sum())
    raw[~identified] = np.nan
    if args.pr_matrix_file:
        check_matches_pr_matrix(dense, identified, args.pr_matrix_file)
    uns = modality_uns(stats, linear=True, quantity_source=DIANN_SOURCES["precursors"].format(q=args.matrix_qvalue))
    uns.update({
        "layers": {"raw": DIANN_SOURCES["precursors.raw"]},
        "n_precursors_outside_pr_matrix": int((~keep).sum()),
        "n_raw_cells_outside_pr_matrix": n_masked,
    })
    return dense, {"raw": raw}, precursors.var.loc[var_names].copy(), uns


def check_matches_pr_matrix(dense: np.ndarray, identified: np.ndarray, path) -> None:
    """G18: same precursor count, same non-empty cells count and the same values as DIA-NN's pr_matrix.

    qpx renames Precursor.Id to peptidoform|charge, so the comparison is by counts and by the
    sorted values (Precursor.Normalised), not by name.
    """
    table = pd.read_csv(path, sep="\t")
    runs = [c for c in table.columns if run_file_stem(c) != c]
    values = table[runs].to_numpy(dtype=np.float64)
    n_cells, n_x = int(np.isfinite(values).sum()), int(identified.sum())
    if len(table) != dense.shape[1] or n_cells != n_x:
        raise GuardError(
            "G18", f"precursors have {dense.shape[1]} vars / {n_x} cells but {Path(path).name} has "
            f"{len(table)} precursors / {n_cells} cells",
        )
    expected = np.sort(values[values > 0]).astype(np.float32)
    got = np.sort(dense[~np.isnan(dense)])
    if expected.shape != got.shape or not np.allclose(got, expected, rtol=1e-6, atol=0):
        raise GuardError("G18", f"precursors.X values differ from {Path(path).name} beyond float32 tolerance")


def assemble_diann(args, mdata, dense):
    """DIA-NN: obs are run stems in both modalities; check them against the SDRF."""
    expected = sdrf_run_stems(args.sdrf_file)
    sdrf = sdrf_obs_table(args.sdrf_file)
    out = {}
    for name, build in (("precursors", diann_precursors), ("proteins", diann_proteins)):
        adata = mdata.mod[name]
        check_obs_names(name, adata.obs_names, expected)
        check_sample_accession(name, adata.obs)
        X, layers, var, uns = build(args, adata, dense[name][0])
        out[name] = {"X": X, "layers": layers, "var": var, "uns": uns,
                     "obs": join_obs_metadata(adata.obs, sdrf, "sdrf_")}
    return out, "run_stem"


def fragpipe_protein_var(pg_path, var: pd.DataFrame) -> pd.DataFrame:
    """Add protein_accession / protein_group / is_contaminant from combined_protein.tsv."""
    table = pd.read_csv(pg_path, sep="\t", dtype=str, keep_default_na=False)
    table.index = table["Protein"].str.split(",").str[0].str.strip()
    rows = table.loc[list(map(str, var.index))]
    extra = rows["Indistinguishable Proteins"] if "Indistinguishable Proteins" in rows.columns else [""] * len(rows)
    groups = [[leader, *split_members(other, ",")] for leader, other in zip(rows.index, extra)]
    return var.join(protein_var_columns(list(rows.index), groups).set_index(var.index))


def assemble_fragpipe(args, mdata, dense):
    """FragPipe: align precursors (keyed by experiment) with proteins (keyed by run)."""
    exp_to_runs = experiment_runs(args)
    protein_source = precheck_fragpipe_columns(read_tsv_header(args.pg_file), exp_to_runs)
    expected_runs = sdrf_run_stems(args.sdrf_file)
    check_obs_names("experiment->run mapping", [r for runs in exp_to_runs.values() for r in runs], expected_runs)
    proteins, precursors = mdata.mod["proteins"], mdata.mod["precursors"]
    check_obs_names("proteins", proteins.obs_names, expected_runs)
    check_obs_names("precursors", precursors.obs_names, exp_to_runs)
    check_sample_accession("proteins", proteins.obs)

    protein_obs = join_obs_metadata(proteins.obs, sdrf_obs_table(args.sdrf_file), "sdrf_")
    if args.experiment_annotation_file:
        protein_obs = join_obs_metadata(protein_obs, annotation_obs_table(args.experiment_annotation_file), "annotation_")
    protein_matrix = dense["proteins"][0]
    fractionated = any(len(runs) > 1 for runs in exp_to_runs.values())
    if fractionated:
        check_fraction_rows_identical(protein_matrix, proteins.obs_names, exp_to_runs)
        protein_matrix, protein_obs = collapse_fragpipe_obs(protein_matrix, protein_obs, exp_to_runs)
    # G8 against combined_protein.tsv: qpx drops source zeros, so count them in the TSV
    source_missing = fragpipe_protein_source_missing(args.pg_file, exp_to_runs, proteins.var_names)
    check_nan_count("proteins", protein_matrix, sum(source_missing.values()))
    protein_stats = {**dense["proteins"][1], **source_missing}
    protein_uns = modality_uns(protein_stats, linear=True, quantity_source=protein_source)
    protein_uns.update(fragpipe_dropped_proteins(args.pg_file, sorted(exp_to_runs), proteins.var_names))
    order, precursor_obs = align_fragpipe_precursor_obs(precursors.obs_names, protein_obs, exp_to_runs, fractionated)
    out = {
        "precursors": {
            "X": dense["precursors"][0][order], "layers": {}, "var": precursors.var.copy(), "obs": precursor_obs,
            "uns": modality_uns(dense["precursors"][1], linear=True, quantity_source="<exp> Intensity (combined_ion.tsv)"),
        },
        "proteins": {
            "X": protein_matrix, "layers": {}, "var": fragpipe_protein_var(args.pg_file, proteins.var.copy()),
            "obs": protein_obs, "uns": protein_uns,
        },
    }
    return out, "experiment" if fractionated else "run_stem"


def build_output(mdata, parts):
    """Final MuData with dense NaN matrices and provenance; runs G4, G5, G9 (final matrices), G10."""
    import anndata as ad
    import mudata as mu

    check_aligned_obs({name: parts[name]["obs"].index for name in MODALITIES})
    modalities = {}
    for name in MODALITIES:
        part = parts[name]
        check_nonempty(name, part["X"])
        check_sample_accession(name, part["obs"])
        check_value_range(name, part["X"], linear=True)
        for layer, matrix in part["layers"].items():
            check_value_range(f"{name}.layers['{layer}']", matrix, linear=True)
        layers = {k: v.astype(np.float32) for k, v in part["layers"].items()}
        modalities[name] = ad.AnnData(X=part["X"].astype(np.float32), obs=part["obs"].copy(), var=part["var"],
                                      uns=part["uns"], layers=layers)
    out = mu.MuData(modalities)
    out.uns.update(dict(mdata.uns))
    if "feature_mapping" in mdata.varp:
        if out.var_names.equals(mdata.var_names):
            out.varp["feature_mapping"] = mdata.varp["feature_mapping"]
        else:  # DIA-NN proteins were restricted to pg_matrix: keep the matching rows / columns
            keep = mdata.var_names.get_indexer(out.var_names)
            if (keep < 0).any():
                raise GuardError("G10", "var axis gained names while rebuilding the MuData; feature_mapping cannot be carried")
            out.varp["feature_mapping"] = mdata.varp["feature_mapping"][keep][:, keep]
    return out


def export(args) -> None:
    if args.engine == "fragpipe":
        run_precheck(args)  # G3 first, so a missing MaxLFQ column is named as the cause
    mdata = build_qpx_mudata(args.qpx_folder)
    check_failed_modalities(mdata.uns)
    hints = MISSING_MODALITY_HINT[args.engine]
    check_modality_present(mdata.mod, "proteins", "G2", hints["proteins"])
    check_modality_present(mdata.mod, "precursors", "G7", hints["precursors"])
    dense = densify_modalities(args, mdata)
    assemble = assemble_diann if args.engine == "diann" else assemble_fragpipe
    parts, obs_key = assemble(args, mdata, dense)
    out = build_output(mdata, parts)

    summary = {"engine": args.engine, "prefix": args.prefix, "obs_key": obs_key, "modalities": ",".join(out.mod)}
    summary["guards_passed"] = GUARDS_RUN[args.engine] + (",G18" if args.pr_matrix_file else "")
    for name in MODALITIES:
        adata = out.mod[name]
        summary.update(summary_rows("", name, adata.X, adata.uns, {k: adata.layers[k] for k in adata.layers.keys() if k is not None}))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out.write(outdir / f"{args.prefix}.h5mu")
    # the standalone .h5ad has no MuData wrapper: carry the dataset-level provenance
    # (qpx writer_version, software_name, qpx_failed_modalities, ...) into its own uns
    proteins = out.mod["proteins"]
    proteins.uns = {**out.uns, **proteins.uns}
    proteins.write_h5ad(outdir / f"{args.prefix}.proteins.h5ad")
    write_summary_tsv(outdir / f"{args.prefix}.summary.tsv", summary)
    print(f"wrote {args.prefix}.h5mu: " + ", ".join(f"{n} {out.mod[n].shape}" for n in out.mod))


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        if args.precheck:
            print(f"precheck OK: protein quantity source {run_precheck(args)}")
        else:
            export(args)
    except GuardError as error:
        print(str(error), file=sys.stderr)
        return GUARD_EXIT_CODE
    return 0


if __name__ == "__main__":
    sys.exit(main())
