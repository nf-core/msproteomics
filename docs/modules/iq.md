# IQ

Protein quantification using the iq R package with MaxLFQ algorithm.

## Description

IQ runs the `iq` R package to perform MaxLFQ protein quantification from DIA-NN report output.
It removes contaminant proteins matching a user-defined regex pattern and performs label-free quantification using an R template script (`iq.r`).
Used in the main DIA workflow for protein-level quantification after DIA-NN analysis.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| report_path | path | DIA-NN main report TSV file (`report.tsv`) |
| q | val | Precursor-level FDR threshold (q-value cutoff) |
| pgq | val | Protein group-level FDR threshold (q-value cutoff) |
| contaminant_pattern | val | Regex pattern to identify contaminant proteins (e.g., `Cont_`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| contaminants_removed_report | path | DIA-NN report with contaminant entries removed (`contaminants_removed_*.tsv`) |
| maxlfq | path | MaxLFQ protein quantification matrix (`maxlfq.tsv`) |
| versions | path | Software versions (`versions.yml`) for r-base and iq |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| (none) | - | Logic is in the R template `iq.r`; no ext.args used |

## Container

`docker.io/dongzehe/iq:0.0.1`

## Usage

Used in the main workflow (`workflows/msproteomics.nf`) as part of the DIA analysis path for MaxLFQ protein quantification from DIA-NN output.

## References

- [iq R package GitHub](https://github.com/tvpham/iq)
