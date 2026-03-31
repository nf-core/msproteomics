# MSSTATS_LFQ

Statistical analysis of label-free quantification data using MSstats.

## Description

MSSTATS_LFQ performs statistical analysis of quantitative mass spectrometry-based proteomics experiments using the MSstats Bioconductor package.
It runs `msstats_plfq.R` with configurable contrast definitions, reference conditions, and quantification parameters to produce CSV results and optional PDF plots.
This module fits at the end of the DIA pipeline for downstream statistical analysis of quantified proteins.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| msstats_csv_input | path | MSstats-formatted CSV input file (`out_msstats.csv`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| msstats_csv | path | MSstats result CSV files (`*.csv`) |
| log | path | Log file (`*.log`) |
| versions | path | Software versions file (`versions.yml`) with r-base and bioconductor-msstats versions |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Passed as last positional argument to `msstats_plfq.R` |
| params.contrasts | - | Contrast definitions for statistical testing |
| params.ref_condition | - | Reference condition for contrasts |
| params.msstats_remove_one_feat_prot | - | Remove single-feature proteins |
| params.msstatslfq_removeFewMeasurements | - | Remove features with few measurements |
| params.msstatslfq_feature_subset_protein | - | Feature subset for protein quantification |
| params.msstatslfq_quant_summary_method | - | Quantification summary method |
| params.msstats_threshold | - | Significance threshold |

## Container

`biocontainers/bioconductor-msstats:4.14.0--r44he5774e6_0`

## Usage

Used in the main `msproteomics.nf` workflow (`workflows/msproteomics.nf`) for statistical downstream analysis of DIA label-free quantification results.
PDF plot outputs are noted as unstable.

## References

- [MSstats GitHub](https://github.com/Vitek-Lab/MSstats)
- [MSstats Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html)
- [MSstats Documentation](https://www.bioconductor.org/packages/release/bioc/manuals/MSstats/man/MSstats.pdf)
