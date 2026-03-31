# TMT_LABELCHECK_ANALYZE

Analyzes TMT labeling efficiency from PSM data with per-peptide site-level calculation.

## Description

This module calculates TMT labeling efficiency based on (Labeled Sites / Total Labelable Sites) where labelable sites = lysine residues + N-terminus per peptide.
It properly handles N-terminal acetylation (42.0106 Da) as biological blocking rather than labeling failure.
Supports two modes: `psm` (finds and concatenates per-sample psm.tsv files, runs `tmt_qc.py analyze`) and `ionquant` (runs `tmt_qc_ionquant.py` on combined_modified_peptide.tsv).
Generates QC reports in HTML, Markdown, and TSV formats with pass/warn/fail thresholds.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| input_files | path | Per-sample result directories (psm mode) or combined_modified_peptide.tsv (ionquant mode) |
| tmt_type | val(string) | TMT reagent type (TMT0, TMT2, TMT6, TMT10, TMT11, TMT16, TMT18, TMT35, TMTPRO) |
| mode | val(string) | `'psm'` for per-sample PSM files, `'ionquant'` for combined_modified_peptide.tsv |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| html_report | tuple(val, path) | `report.html` rich HTML report with color-coded labeling status |
| md_report | tuple(val, path) | `report.md` Markdown report for documentation |
| summary | tuple(val, path) | `labeling_summary.tsv` machine-readable summary |
| per_sample | tuple(val, path) | `per_sample_efficiency.csv` per-sample efficiency data |
| versions_python | topic: versions | Python version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments for tmt_qc scripts |
| ext.pass_threshold | 0.95 | Labeling efficiency pass threshold |
| ext.warn_threshold | 0.85 | Labeling efficiency warning threshold |
| ext.prob_threshold | 0.95 | PSM probability threshold for filtering |

## Container

`community.wave.seqera.io/library/matplotlib_pandas:76f3c63ec67531f0`

## Usage

Used in the `TMT_LABELCHECK` subworkflow (`subworkflows/local/tmt_labelcheck/main.nf`) as the final QC step after database search and filtering, to verify TMT labeling efficiency before proceeding with quantification.

## References

- [FragPipe](https://fragpipe.nesvilab.org/)
