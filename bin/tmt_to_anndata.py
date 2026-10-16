#!/usr/bin/env python3
"""Convert TMT-Integrator abundance reports to guarded .h5mu / .h5ad files.

Input: every ``abundance_<level>_<norm>.tsv`` TMT-Integrator wrote (levels gene,
protein, peptide, single-site, multi-site; norms None, MD, GN, SL+IRS) and the
TMT annotation (``plex channel sample sample_name condition replicate``).
Header layouts follow TMT-Integrator ``ReportGenerator.writeTitle`` /
``ReportInfo`` (data/TMT-Integrator/src/main/java/tmtintegrator/):
annotation columns, ``ReferenceIntensity``, one column per sample, then
optional ``RefInt_*`` / ``RefDInt_*`` columns.

Output per normalization ``<norm>``:

    <prefix>.tmt_<norm>.h5mu                one modality per level (gene, protein, peptide, single_site, multi_site)
    <prefix>.tmt_<level>_<norm>.h5ad        one per level
    <prefix>.tmt.summary.tsv                sorted key<TAB>value rows over every modality

obs = samples (annotation ``sample_name``), var = the report ``Index`` (numeric
header columns as numbers, plus ``protein_accession``, ``protein_group`` and
``is_contaminant``), X = dense float32. Missing values are NaN: ``NA``/empty
always; on a linear scale also values <= 0. On log2 scale 0 is kept.

Scale per normalization (TMT-Integrator source, data/TMT-Integrator/src/main/java/tmtintegrator/):
None / MD / GN follow ``log2transformed`` (``ReportGenerator.writeRatioAbnValues``);
SL+IRS is log2-like whatever the flag: ``PsmNormalizer.ratioToAbundance``
(lines 120-141) converts to ``log2(2^ratio * avgAbundance)``, SL/IRS then scale
those values, and ``ReportGenerator.writeValues`` (344-362) writes them through
``writeAbnValues`` without the flag. Its ``ReferenceIntensity`` is linear
(``writeRefInt``, 316-325).

Guards (all exit 3 with ``ERROR: G<n>: ...`` before anything is written):
G4 non-empty, G8 NaN count equals the source missing count, G9 no inf / no
value <= 0 on linear scale (final X), G12 header layout, numeric columns and
sample columns map to the annotation, G13 obs identical across levels of one
normalization, G14 at least one recognised abundance file, G15 the values and
``ReferenceIntensity`` agree with the scale rule (both directions).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from qpx_export_lib import (
    GUARD_EXIT_CODE,
    GuardError,
    check_nan_count,
    check_nonempty,
    check_value_range,
    densify_with_nan,
    modality_uns,
    protein_var_columns,
    split_members,
    summary_rows,
    write_summary_tsv,
)

# ReportInfo.HEADER / HEADER_PROTEIN and ReportGenerator.writeTitle (GENE) in TMT-Integrator.
HEADER_SITE = [
    "Index", "Gene", "ProteinID", "Peptide", "SequenceWindow", "Start", "End", "MaxPepProb",
    "Spectrum Number", "Protein", "Entry Name", "Protein Description", "Mapped Genes", "Mapped Proteins",
]
LEVEL_HEADERS = {
    "gene": ["Index", "NumberPSM", "ProteinID", "MaxPepProb"],
    "protein": [
        "Index", "NumberPSM", "Gene", "MaxPepProb", "Protein", "Protein ID", "Entry Name",
        "Protein Description", "Organism", "Indistinguishable Proteins",
    ],
    "peptide": HEADER_SITE,
    "single-site": HEADER_SITE,
    "multi-site": HEADER_SITE,
}
MODALITY_NAME = {"gene": "gene", "protein": "protein", "peptide": "peptide",
                 "single-site": "single_site", "multi-site": "multi_site"}
FILE_PATTERN = re.compile(r"^abundance_(?P<level>.+)_(?P<norm>None|MD|GN|SL\+IRS)\.tsv$")
REFERENCE_COLUMN = "ReferenceIntensity"
REFINT_PREFIXES = ("RefInt_", "RefDInt_")
MISSING_TOKENS = ("NA", "")
# A log2 abundance above this is 2^64 on the linear scale: the table is not log2.
# Conversely, linear TMT abundances are reporter intensities (~1e3-1e9); a linear
# table whose largest value is below this is a log2 table.
LOG2_MAX_PLAUSIBLE = 64.0
NUMERIC_VAR_COLUMNS = ("NumberPSM", "MaxPepProb", "ReferenceIntensity", "Start", "End", "Spectrum Number")
# (leader column, comma-separated member column) for protein_accession / protein_group per level
PROTEIN_COLUMNS = {
    "gene": ("ProteinID", None),
    "protein": ("Protein", "Indistinguishable Proteins"),
    "peptide": ("Protein", "Mapped Proteins"),
    "single-site": ("Protein", "Mapped Proteins"),
    "multi-site": ("Protein", "Mapped Proteins"),
}
SL_IRS = "SL+IRS"
ANNOTATION_COLUMNS = ["plex", "channel", "sample", "sample_name", "condition", "replicate"]


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--abundance", nargs="+", required=True, help="abundance_<level>_<norm>.tsv files")
    parser.add_argument("--annotation", required=True, help="TMT annotation TSV (plex channel sample sample_name ...)")
    parser.add_argument("--log2transformed", required=True, choices=["true", "false"],
                        help="TMT-Integrator log2transformed setting the reports were written with")
    parser.add_argument("--prefix", required=True, help="output file prefix")
    parser.add_argument("--outdir", default=".", help="output directory (default: .)")
    return parser.parse_args(argv)


def classify_files(paths) -> dict[str, dict[str, Path]]:
    """G14: group abundance files as ``{norm: {level: path}}``; unknown names fail."""
    if not paths:
        raise GuardError("G14", "TMT-Integrator produced no abundance_*.tsv files")
    groups: dict[str, dict[str, Path]] = {}
    for raw in paths:
        path = Path(raw)
        match = FILE_PATTERN.match(path.name)
        if not path.is_file() or match is None or match["level"] not in LEVEL_HEADERS:
            raise GuardError(
                "G14",
                f"'{raw}' is not a readable abundance_<level>_<norm>.tsv with level in {sorted(LEVEL_HEADERS)} "
                "and norm in None/MD/GN/SL+IRS",
            )
        groups.setdefault(match["norm"], {})[match["level"]] = path
    return groups


def read_annotation(path) -> pd.DataFrame:
    """Annotation rows with a sample_name, indexed by sample_name (G12 on bad layout)."""
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    table.columns = [c.strip() for c in table.columns]
    missing = [c for c in ANNOTATION_COLUMNS if c not in table.columns]
    if missing:
        raise GuardError("G12", f"annotation {Path(path).name} lacks column(s) {missing}")
    table = table[ANNOTATION_COLUMNS].apply(lambda col: col.str.strip())
    table = table[table["sample_name"] != ""]
    duplicated = sorted(table.loc[table["sample_name"].duplicated(), "sample_name"].unique())
    if duplicated:
        raise GuardError("G12", f"annotation sample_name values are not unique: {duplicated}")
    return table.set_index("sample_name")


def split_columns(path: Path, level: str, columns: list[str]) -> tuple[list[str], list[str], list[str]]:
    """G12: check the header layout; return (annotation, sample, RefInt) columns."""
    layout = LEVEL_HEADERS[level]
    n = len(layout)
    if columns[:n] != layout or len(columns) <= n or columns[n] != REFERENCE_COLUMN:
        raise GuardError(
            "G12",
            f"{path.name}: header {columns[: n + 1]} does not match the TMT-Integrator {level} layout "
            f"{layout + [REFERENCE_COLUMN]}",
        )
    tail = columns[n + 1:]
    samples = [c for c in tail if not c.startswith(REFINT_PREFIXES)]
    refint = [c for c in tail if c.startswith(REFINT_PREFIXES)]
    if not samples:
        raise GuardError("G12", f"{path.name}: no sample columns after {REFERENCE_COLUMN}")
    return layout + [REFERENCE_COLUMN], samples, refint


def check_samples_annotated(path: Path, samples: list[str], annotation: pd.DataFrame) -> None:
    """G12: every sample column maps to exactly one annotation sample_name.

    The reverse is not required: TMT-Integrator omits reference channels.
    """
    duplicated = sorted({s for s in samples if samples.count(s) > 1})
    if duplicated:
        raise GuardError("G12", f"{path.name}: duplicated sample columns {duplicated}")
    unknown = [s for s in samples if s not in annotation.index]
    if unknown:
        raise GuardError(
            "G12", f"{path.name}: sample columns {unknown} are not annotation sample_name values "
            f"{sorted(annotation.index)}",
        )


def parse_values(path: Path, raw: pd.DataFrame, linear: bool) -> tuple[np.ndarray, int]:
    """Numeric matrix with NaN for NA/empty, and the source missing count for G8."""
    stripped = raw.apply(lambda col: col.str.strip())
    missing = stripped.isin(MISSING_TOKENS)
    numeric = stripped.mask(missing).apply(pd.to_numeric, errors="coerce")
    unparsable = int(numeric.isna().to_numpy().sum() - missing.to_numpy().sum())
    if unparsable:
        raise GuardError("G12", f"{path.name}: {unparsable} sample cells are neither numeric nor NA/empty")
    values = numeric.to_numpy(dtype=np.float64)
    expected_missing = int(missing.to_numpy().sum())
    if linear:
        expected_missing += int((values <= 0).sum())
    return values, expected_missing


def norm_scale(norm: str, log2transformed: bool) -> tuple[bool, str]:
    """(values are log2, rule text) for one normalization; see the module docstring."""
    if norm == SL_IRS:
        return True, "SL+IRS: log2-like regardless of log2transformed (PsmNormalizer.java:120-141, ReportGenerator.java:344-362)"
    flag = "true" if log2transformed else "false"
    return log2transformed, f"{norm}: follows log2transformed={flag} (ReportGenerator.writeRatioAbnValues)"


def check_scale(path: Path, what: str, values: np.ndarray, log2: bool, rule: str) -> None:
    """G15: the values agree with the scale rule, in both directions."""
    finite = values[np.isfinite(values)]  # inf is G9's job
    if finite.size == 0:
        return
    peak = float(np.max(np.abs(finite)))
    if log2 and peak > LOG2_MAX_PLAUSIBLE:
        raise GuardError("G15", f"{path.name}: {what} expected log2 ({rule}) but values reach {peak:g}; the report is linear")
    if not log2 and peak <= LOG2_MAX_PLAUSIBLE:
        raise GuardError(
            "G15", f"{path.name}: {what} expected linear ({rule}) but the largest value is {peak:g}; the report is log2"
        )


def numeric_var(path: Path, table: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """G12: known numeric header columns (and RefInt_*) as float; NA/empty -> NaN."""
    out = table[columns].copy()
    for column in [c for c in columns if c in NUMERIC_VAR_COLUMNS or c.startswith(REFINT_PREFIXES)]:
        text = out[column].str.strip()
        missing = text.isin(MISSING_TOKENS)
        values = pd.to_numeric(text.mask(missing), errors="coerce")
        bad = int(values.isna().sum() - missing.sum())
        if bad:
            raise GuardError("G12", f"{path.name}: column '{column}' has {bad} values that are neither numeric nor NA/empty")
        out[column] = values.astype(np.float64)
    return out


def protein_ids(level: str, table: pd.DataFrame) -> pd.DataFrame:
    leader_col, members_col = PROTEIN_COLUMNS[level]
    leaders = table[leader_col].str.strip().tolist()
    others = table[members_col].tolist() if members_col else [""] * len(table)
    groups = [[leader, *split_members(other, ",")] for leader, other in zip(leaders, others)]
    return protein_var_columns(leaders, groups)


def build_level(path: Path, level: str, norm: str, annotation: pd.DataFrame, log2transformed: bool):
    """One AnnData for one abundance report, after all per-file guards."""
    import anndata as ad

    log2, rule = norm_scale(norm, log2transformed)
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    var_cols, samples, refint = split_columns(path, level, list(table.columns))
    check_samples_annotated(path, samples, annotation)
    duplicated = sorted(table.loc[table["Index"].duplicated(), "Index"].unique())
    if duplicated:
        raise GuardError("G12", f"{path.name}: duplicated Index values {duplicated[:5]}")
    values, expected_missing = parse_values(path, table[samples], linear=not log2)
    dense, stats = densify_with_nan(values.T, linear=not log2)
    name = MODALITY_NAME[level]
    check_nan_count(name, dense, expected_missing)
    check_scale(path, "sample columns", dense, log2, rule)
    var = numeric_var(path, table, var_cols[1:] + refint)
    reference_log2 = log2 and norm != SL_IRS  # writeRefInt: SL+IRS reference is linear
    check_scale(path, REFERENCE_COLUMN, var[REFERENCE_COLUMN].to_numpy(), reference_log2, rule)
    check_nonempty(name, dense)
    var = pd.concat([var, protein_ids(level, table)], axis=1).set_index(table["Index"].rename(None))
    obs = annotation.loc[samples].rename_axis("sample_name")
    uns = modality_uns(stats, linear=not log2, quantity_source=path.name)
    uns.update({"log2transformed": log2, "log2transformed_setting": log2transformed, "scale": "log2" if log2 else "linear",
                "scale_rule": rule, "reference_intensity_scale": "log2" if reference_log2 else "linear",
                "prot_norm": norm, "groupby": level})
    return ad.AnnData(X=dense, obs=obs, var=var, uns=uns)


def check_same_obs(norm: str, modalities: dict) -> None:
    """G13: every level of one normalization has the same samples in the same order."""
    names = list(modalities)
    ref = list(modalities[names[0]].obs_names)
    for name in names[1:]:
        if list(modalities[name].obs_names) != ref:
            raise GuardError(
                "G13",
                f"norm {norm}: '{name}' samples {list(modalities[name].obs_names)} differ from '{names[0]}' samples {ref}",
            )


def convert(args) -> None:
    import mudata as mu

    log2transformed = args.log2transformed == "true"
    groups = classify_files(args.abundance)
    annotation = read_annotation(args.annotation)
    outputs = {}
    summary = {"engine": "tmt", "prefix": args.prefix, "log2transformed": args.log2transformed,
               "norms": ",".join(sorted(groups)), "guards_passed": "G4,G8,G9,G12,G13,G14,G15"}
    for norm in sorted(groups):
        levels = [lvl for lvl in LEVEL_HEADERS if lvl in groups[norm]]
        modalities = {
            MODALITY_NAME[lvl]: build_level(groups[norm][lvl], lvl, norm, annotation, log2transformed) for lvl in levels
        }
        check_same_obs(norm, modalities)
        for name, adata in modalities.items():  # G9 on the final matrices
            check_value_range(f"tmt_{norm}.{name}", adata.X, linear=adata.uns["scale"] == "linear")
        outputs[norm] = modalities
        summary[f"tmt_{norm}.modalities"] = ",".join(modalities)
        for name, adata in modalities.items():
            summary.update(summary_rows(f"tmt_{norm}.", name, adata.X, adata.uns))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for norm, modalities in outputs.items():
        mdata = mu.MuData(modalities)
        scale = modalities[next(iter(modalities))].uns["scale"]
        mdata.uns.update({"log2transformed": scale == "log2", "log2transformed_setting": log2transformed,
                          "scale": scale, "prot_norm": norm, "source": "TMT-Integrator"})
        mdata.write(outdir / f"{args.prefix}.tmt_{norm}.h5mu")
        for name, adata in modalities.items():  # standalone .h5ad: carry the dataset-level uns
            adata.uns = {**mdata.uns, **adata.uns}
            adata.write_h5ad(outdir / f"{args.prefix}.tmt_{name}_{norm}.h5ad")
    write_summary_tsv(outdir / f"{args.prefix}.tmt.summary.tsv", summary)
    for norm, modalities in outputs.items():
        print(f"wrote {args.prefix}.tmt_{norm}.h5mu: " + ", ".join(f"{n} {a.shape}" for n, a in modalities.items()))


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        convert(args)
    except GuardError as error:
        print(str(error), file=sys.stderr)
        return GUARD_EXIT_CODE
    return 0


if __name__ == "__main__":
    sys.exit(main())
