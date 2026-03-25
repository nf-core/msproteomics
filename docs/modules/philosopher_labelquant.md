# PHILOSOPHER_LABELQUANT

Extracts TMT/isobaric reporter ion intensities from mzML spectra using Philosopher labelquant.

## Description

This module runs Philosopher `labelquant` to extract TMT/isobaric reporter ion intensities from mzML spectra.
It runs inside the results directory workspace (containing `.meta/` from PHILOSOPHER_FILTER or FREEQUANT), reading `.meta/` binaries and mzML files.
The updated `.meta/` workspace is passed to downstream PHILOSOPHER_REPORT, which regenerates TSVs with TMT intensity columns.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| results_dir | path | Directory with psm.tsv and `.meta/` workspace from PHILOSOPHER_FILTER or FREEQUANT |
| mzml_dir | path | Directory containing mzML files for the experiment/group |
| annot | path | TMT annotation file mapping channels to samples |
| config_cli | val(string) | Labelquant CLI flags (`--tol`, `--level`, `--plex`, `--brand`, `--minprob`, `--purity`, `--removelow`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Output directory `${prefix}` with updated `.meta/` containing TMT intensities |
| versions_philosopher | topic: versions | Philosopher version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments appended after `config_cli` and `--annot`/`--dir` flags |
| ext.prefix | `${meta.id}` | Output directory name |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_QUANT` subworkflow (`subworkflows/local/fragpipe_quant/main.nf`) for TMT reporter ion extraction.
The output `.meta/` workspace is then passed to PHILOSOPHER_REPORT to generate final TSV files with quantification data.

## References

- [Philosopher](https://philosopher.nesvilab.org/)
- [Philosopher wiki](https://github.com/Nesvilab/philosopher/wiki)
- [Publication](https://doi.org/10.1038/s41592-020-0912-y)
