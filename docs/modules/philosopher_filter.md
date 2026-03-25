# PHILOSOPHER_FILTER

Performs FDR filtering and report generation for proteomics results using Philosopher.

## Description

This module runs Philosopher `filter` and optionally `report` for FDR filtering at PSM, peptide, and protein levels.
It annotates the database, runs filter with specified flags (e.g., `--sequential --razor --prot 0.01 --picked`), then optionally runs report to generate TSV files.
When `skip_report=true`, the `.meta/` workspace is preserved for downstream tools (FreeQuant, LabelQuant) and placeholder TSVs are created.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| pepxml_files | path | PepXML files from search/validation (`*.pep.xml`) |
| protxml | path | ProtXML from ProteinProphet (`*.prot.xml`) |
| filter_flags | val(string) | Filter CLI flags (e.g., `--sequential --razor --prot 0.01 --picked`) |
| report_flags | val(string) | Report CLI flags (e.g., `--msstats`) |
| fasta | path | Philosopher database with decoys/contaminants from PHILOSOPHER_DATABASE |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Output directory `${prefix}` containing all result files |
| psms | tuple(val, path) | `${prefix}/psm.tsv` PSM-level results |
| peptides | tuple(val, path) | `${prefix}/peptide.tsv` peptide-level results |
| proteins | tuple(val, path) | `${prefix}/protein.tsv` protein-level results |
| ions | tuple(val, path) | `${prefix}/ion.tsv` ion-level results (optional) |
| versions_philosopher | topic: versions | Philosopher version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments for `philosopher filter` |
| ext.prefix | `${meta.id}` | Output directory name |
| ext.decoy_tag | `rev_` | Decoy protein prefix |
| ext.skip_report | false | When true, preserve `.meta/` workspace for downstream FreeQuant/LabelQuant |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_INFERENCE` subworkflow (`subworkflows/local/fragpipe_inference/main.nf`) for per-sample FDR filtering after ProteinProphet protein inference.
In TMT workflows with `skip_report=true`, the `.meta/` workspace is passed to PHILOSOPHER_LABELQUANT and then PHILOSOPHER_REPORT.

## References

- [Philosopher](https://philosopher.nesvilab.org/)
- [Philosopher wiki](https://github.com/Nesvilab/philosopher/wiki)
- [Publication](https://doi.org/10.1038/s41592-020-0912-y)
