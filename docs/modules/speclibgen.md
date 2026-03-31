# SPECLIBGEN

Generates spectral libraries using EasyPQP via FragPipe-SpecLib.

## Description

SpecLibGen generates spectral libraries from database search results using EasyPQP (Easy Peptide Quantification Pipeline) within the FragPipe framework.
It performs per-sample PSM conversion (`convertpsm`), then builds a unified spectral library (`library`), supporting RT calibration (noiRT, ciRT, Pierce_iRT, Biognosys_iRT, or custom file) and IM calibration.
Output formats include TSV, Spectronaut (`.speclib`), and Parquet, with automatic fallback to ciRT if RT alignment fails.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| psm_files | path (stageAs `psms/sample_*`) | Per-sample psm.tsv files from PHILOSOPHER_FILTER |
| peptide_files | path (stageAs `peps/sample_*`) | Per-sample peptide.tsv files from PHILOSOPHER_FILTER |
| mzml_files | path (stageAs `spectra/*`) | Mass spec files for spectrum extraction |
| config_cli | val(string) | Bash-sourceable KEY='value' lines (CONVERT_ARGS, FRAGMENT_TYPES, LIBRARY_ARGS, RT_CAL, IM_CAL, KEEP_INTERMEDIATE) |
| fasta | path | Protein sequence database |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Output directory `${prefix}` containing all spectral library files |
| library_tsv | tuple(val, path) | `${prefix}/library.tsv` spectral library in TSV format (optional) |
| library_speclib | tuple(val, path) | `${prefix}/library.speclib` in speclib binary format (optional) |
| library_parquet | tuple(val, path) | `${prefix}/library.parquet` in Parquet columnar format (optional) |
| versions_speclibgen | topic: versions | SpecLibGen version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments appended to EasyPQP commands |
| ext.prefix | `${meta.id}` | Output directory name |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory for binary discovery |
| ext.decoy_tag | `rev_` | Decoy protein prefix |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_SPECLIB` subworkflow (`subworkflows/local/fragpipe_speclib/main.nf`) for generating spectral libraries from DDA search results, which can then be used for DIA data analysis or targeted proteomics.

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
- [SpecLib tutorial](https://fragpipe.nesvilab.org/docs/tutorial_speclib.html)
