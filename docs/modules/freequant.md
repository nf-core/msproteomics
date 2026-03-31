# FREEQUANT

Perform label-free quantification using Philosopher's `freequant` command.

## Description

FreeQuant extracts ion intensity information from mass spectrometry data files to quantify peptide and protein abundances.
It initializes a Philosopher workspace in the results directory, runs FreeQuant with the provided CLI config, and preserves the `.meta/` workspace for downstream PHILOSOPHER_REPORT.
Used within the FRAGPIPE_QUANT subworkflow as an alternative to IonQuant for label-free quantification.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| results_dir | path | Directory containing `psm.tsv` from PHILOSOPHER_FILTER |
| mzml_dir | path | Directory containing mass spectrometry files (mzML or RAW) |
| config_cli | val(string) | FreeQuant parameters as CLI string (e.g., `--ptw 0.4 --tol 10`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(meta, path) | Output directory with results and `.meta/` workspace |
| ions | tuple(meta, path) | Ion-level quantification results (`ion.tsv`, optional) |
| versions_philosopher | tuple | Philosopher software version (topic channel) |
| versions_fragpipe | tuple | FragPipe software version (topic channel) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the `philosopher freequant` command |
| ext.prefix | meta.id | Override for the output directory name |
| ext.fragpipe_tools_dir | (auto-detected) | Override for the FragPipe tools directory path |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_QUANT` subworkflow (`subworkflows/local/fragpipe_quant/main.nf`) for label-free ion intensity quantification.
Automatically detects `.raw` files and adds the `--raw` flag.

## References

- [Philosopher homepage](https://philosopher.nesvilab.org/)
- [Philosopher GitHub](https://github.com/Nesvilab/philosopher)
- [Philosopher Wiki](https://github.com/Nesvilab/philosopher/wiki)
