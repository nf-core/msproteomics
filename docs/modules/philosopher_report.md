# PHILOSOPHER_REPORT

Generates TSV reports from Philosopher workspace binaries.

## Description

This module runs Philosopher `report` to generate TSV result files (psm.tsv, peptide.tsv, protein.tsv, ion.tsv) from `.meta/` workspace binaries.
It is used after PHILOSOPHER_FILTER (with `skip_report=true`) followed by FreeQuant or LabelQuant, which update the `.meta/` binaries with quantification data.
The module re-initializes the workspace and re-annotates the database to ensure accessibility, since the workspace may come from a different task's work directory.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| results_dir | path | Directory containing `.meta/` workspace binaries from PHILOSOPHER_FILTER/FREEQUANT/LABELQUANT |
| report_cli | val(string) | Report CLI flags (e.g., `--msstats`) |
| fasta | path | Philosopher database with decoys/contaminants from PHILOSOPHER_DATABASE |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Output directory `${prefix}` containing all TSV files |
| psms | tuple(val, path) | `${prefix}/psm.tsv` PSM-level results with quantification intensities |
| peptides | tuple(val, path) | `${prefix}/peptide.tsv` peptide-level results |
| proteins | tuple(val, path) | `${prefix}/protein.tsv` protein-level results |
| ions | tuple(val, path) | `${prefix}/ion.tsv` ion-level results (optional) |
| versions_philosopher | topic: versions | Philosopher version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments for `philosopher report` |
| ext.prefix | `${meta.id}` | Output directory name |
| ext.decoy_tag | `rev_` | Decoy protein prefix |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_QUANT` subworkflow (`subworkflows/local/fragpipe_quant/main.nf`) after PHILOSOPHER_LABELQUANT to generate final TSV files with correct TMT reporter ion intensities.

## References

- [Philosopher](https://philosopher.nesvilab.org/)
- [Philosopher wiki](https://github.com/Nesvilab/philosopher/wiki)
- [Publication](https://doi.org/10.1038/s41592-020-0912-y)
