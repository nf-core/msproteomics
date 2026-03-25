# DIAUMPIRE

Generate pseudo-spectra from DIA mass spectrometry data.

## Description

DIA-Umpire Signal Extraction (SE) deconvolves multiplexed DIA spectra to generate pseudo-MS/MS spectra that can be searched with conventional database search engines.
It produces three quality tiers of pseudo-spectra (Q1, Q2, Q3) in mzML and/or MGF format, and falls back to default parameters if no config file is provided.
Used within the FRAGPIPE_CONVERT subworkflow for DIA data preprocessing.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| dia_file | path | DIA mass spectrometry data file (`*.mzML`, `*.mzXML`, or `*.raw`) |
| config_file | path | DIA-Umpire parameters file (`*.params`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(meta, path) | Directory containing DIA-Umpire results and log |
| q1_spectra | tuple(meta, path) | Highest quality tier pseudo-spectra (`*_Q1.mzML`, optional) |
| q2_spectra | tuple(meta, path) | Medium quality tier pseudo-spectra (`*_Q2.mzML`, optional) |
| q3_spectra | tuple(meta, path) | Lowest quality tier pseudo-spectra (`*_Q3.mzML`, optional) |
| mgf_spectra | tuple(meta, path) | Pseudo-spectra in MGF format (`*_Q*.mgf`, optional) |
| versions_diaumpire | tuple | DIA-Umpire software version (topic channel) |
| versions_fragpipe | tuple | FragPipe software version (topic channel) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the DIA-Umpire SE Java command |
| ext.prefix | meta.id | Override for the output directory name |
| ext.fragpipe_tools_dir | (auto-detected) | Override for the FragPipe tools directory path |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_CONVERT` subworkflow (`subworkflows/local/fragpipe_convert/main.nf`) to generate pseudo-spectra from DIA data.
Enabled conditionally via `shouldRunTool(configs, 'diaumpire')`.

## References

- [DIA-Umpire homepage](https://diaumpire.nesvilab.org/)
- [DIA-Umpire GitHub](https://github.com/Nesvilab/DIA-Umpire)
- [Publication: doi:10.1038/nmeth.3255](https://doi.org/10.1038/nmeth.3255)
