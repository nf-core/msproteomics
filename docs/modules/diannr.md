# DIANNR

R-based post-processing of DIA-NN output for contaminant removal and MaxLFQ quantification.

## Description

This module uses the `diann-rpackage` for post-processing DIA-NN report TSV files.
It removes contaminant proteins, applies q-value and protein-group q-value filters, and generates precursor/peptide/protein-level quantification tables with MaxLFQ normalization.
The module uses an R template script (`diannr.r`) that receives parameters via Nextflow input variables.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| report_path | path | DIA-NN main report TSV file |
| q | val | Precursor-level FDR threshold (q-value cutoff) |
| pgq | val | Protein group-level FDR threshold (q-value cutoff) |
| contaminant_pattern | val | Regex pattern to identify contaminant proteins (e.g., `"Cont_"`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| contaminants_removed_report | path | `contaminants_removed_*.tsv` DIA-NN report with contaminants removed |
| precursors | path | `precursors.tsv` precursor-level quantification matrix |
| peptides | path | `peptides.tsv` peptide-level quantification matrix |
| peptides_maxlfq | path | `peptides_maxlfq.tsv` peptide-level MaxLFQ quantification matrix |
| unique_genes | path | `unique_genes.tsv` gene-level unique peptide quantification matrix |
| protein_groups_maxlfq | path | `protein_groups_maxlfq.tsv` protein group-level MaxLFQ quantification matrix |
| versions | path | `versions.yml` with R-base and diann-rpackage versions |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| params.diann_maxlfq_q | - | Q-value threshold passed as the `q` input |
| params.diann_maxlfq_pgq | - | Protein group q-value threshold passed as the `pgq` input |
| params.contaminant_pattern | - | Contaminant pattern passed as the `contaminant_pattern` input |

## Container

`docker.io/dongzehe/diann-r:0.0.1`

## Usage

Used in the main DIA workflow (`workflows/msproteomics.nf`) for post-processing DIA-NN output, generating MaxLFQ-normalized quantification tables at multiple levels (precursor, peptide, protein group).

## References

- [diann-rpackage](https://github.com/vdemichev/diann-rpackage)
- [DIA-NN](https://github.com/vdemichev/DiaNN)
