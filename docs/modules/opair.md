# OPAIR

O-glycoproteomics analysis for O-glycan identification and site localization.

## Description

O-Pair performs O-glycoproteomics analysis by identifying and localizing O-linked glycan modifications on peptides.
It uses paired fragmentation spectra to determine glycan compositions and attachment sites, generating glycoform-level results.
O-Pair is part of the FragPipe glycoproteomics suite and runs as a .NET application (`dotnet CMD.dll`).

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| psm_file | path | PSM file from PHILOSOPHER_FILTER (`psm.tsv`) |
| mzml_files | path | Mass spectrometry data files (`*.mzML`, staged in `spectra/` directory) |
| config_cli | val(string) | O-Pair parameters as CLI string (e.g., `-b 20 -c 20 -n 3`) |
| glycan_db | path | Glycan database file for O-glycan matching |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Results directory containing O-Pair output and logs |
| results | tuple(val, path) | O-Pair glycan identification results (`*_opair_results.tsv`, optional) |
| glycoforms | tuple(val, path) | Glycoform-level results (`*_opair_glycoforms.tsv`, optional) |
| versions_opair | tuple (topic: versions) | O-Pair version |
| versions_fragpipe | tuple (topic: versions) | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to `dotnet CMD.dll` after positional arguments |
| ext.fragpipe_tools_dir | (default path) | Override FragPipe tools directory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_GLYCO` subworkflow (`subworkflows/local/fragpipe_glyco/main.nf`) for O-glycoproteomics analysis.
The `config_cli` parameter provides primary O-Pair parameters (tolerance, glycan settings, oxonium filtering, glycan sites) as a CLI string.

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
- [FragPipe Glycoproteomics Tutorial](https://fragpipe.nesvilab.org/docs/tutorial_glyco.html)
