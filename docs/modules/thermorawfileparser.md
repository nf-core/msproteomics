# THERMORAWFILEPARSER

Converts Thermo `.raw` files to open mass spectrometry formats.

## Description

ThermoRawFileParser is a wrapper around the .NET ThermoFisher ThermoRawFileReader library for running on Linux with Mono.
It converts Thermo `.raw` files to open formats including mzML, mgf, or parquet.
The output format is determined by `ext.args` flags (`--format 0`=mgf, `1`/`2`=mzML, `3`=parquet), and gzip compression is supported via `--gzip`.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata |
| raw | path | Thermo `.raw` file |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| spectra | tuple(val, path) | Converted spectra file (`*.mzML`, `*.mzML.gz`, `*.mgf`, `*.mgf.gz`, `*.parquet`, `*.parquet.gz`) |
| versions | path | `versions.yml` with software versions |
| versions_thermorawfileparser | topic: versions | ThermoRawFileParser version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Command-line arguments for ThermoRawFileParser (controls format and compression) |
| ext.prefix | `${meta.id}` | Output file prefix |

## Container

`quay.io/biocontainers/thermorawfileparser:1.4.5--h05cac1d_1`

## Usage

Used in the `FILE_PREPARATION` subworkflow (`subworkflows/local/file_preparation/main.nf`) and the `FRAGPIPE_CONVERT` subworkflow (`subworkflows/local/fragpipe_convert/main.nf`) for converting Thermo `.raw` files to mzML before database search.
Also used as `THERMORAWFILEPARSER_HEADLESS` in `workflows/msproteomics.nf` for headless raw file conversion.

## References

- [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)
- [Publication](https://doi.org/10.1021/acs.jproteome.9b00328)
