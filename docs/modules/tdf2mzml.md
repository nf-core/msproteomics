# TDF2MZML

Converts Bruker `.d` (TDF) files to mzML format.

## Description

This module converts Bruker `.d` raw files to indexed mzML format using `tdf2mzml.py`.
It renames the output to use a clean base name derived from the input file.
The original `.d` directory is preserved and emitted as an output alongside the converted mzML file.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata (uses `meta.mzml_id` for tag) |
| rawfile | path | Bruker `.d` raw file |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| mzmls_converted | tuple(val, path) | `*.mzML` converted file |
| dotd_files | tuple(val, path) | `*.d` original file (preserved and renamed) |
| versions | path | `versions.yml` with software versions |
| log | path | `*.log` conversion log |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments appended to `tdf2mzml.py -i` |

## Container

`quay.io/bigbio/tdf2mzml@sha256:98a38d2e9d7803e75dfdfdd7282567fec6a4f57df6859c67e0307483c9a710e9`

## Usage

Used in the `FILE_PREPARATION` subworkflow (`subworkflows/local/file_preparation/main.nf`) to convert Bruker `.d` files to mzML when `params.convert_dotd` is enabled.

## References

- [tdf2mzml](https://github.com/mafreitas/tdf2mzml)
