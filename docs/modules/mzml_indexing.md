# MZML_INDEXING

Indexes mzML files using OpenMS FileConverter for efficient random access.

## Description

MZML_INDEXING converts input mzML files to indexed mzML format using OpenMS FileConverter.
Indexed mzML files support efficient random access by downstream tools, which is required for many proteomics analysis steps.
This module is part of the file preparation stage in the pipeline.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map (uses `meta.mzml_id` for tagging) |
| mzmlfile | path | Input mzML file (`*.mzML`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| mzmls_indexed | tuple(val, path) | Indexed mzML file in `out/` directory |
| versions | path | FileConverter version (`versions.yml`) |
| log | path | Log file (`*.log`) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to `FileConverter` command |

## Container

`ghcr.io/bigbio/openms-tools-thirdparty:2025.04.14`

## Usage

Used in the `FILE_PREPARATION` subworkflow (`subworkflows/local/file_preparation/main.nf`) for converting mzML files to indexed format before downstream analysis.

## References

- [OpenMS FileConverter](http://www.openms.de/doxygen/nightly/html/TOPP_FileConverter.html)
- [OpenMS homepage](https://www.openms.de/)
