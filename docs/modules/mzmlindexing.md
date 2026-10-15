# MZMLINDEXING

Indexes mzML files using OpenMS FileConverter for efficient random access.

## Description

MZMLINDEXING converts input mzML files to indexed mzML format using OpenMS FileConverter.
Indexed mzML files support efficient random access by downstream tools, which is required for many proteomics analysis steps.
This module is part of the file preparation stage in the pipeline.

## Inputs

| Channel  | Type     | Description                                           |
| -------- | -------- | ----------------------------------------------------- |
| meta     | val(map) | Sample metadata map (uses `meta.mzml_id` for tagging) |
| mzmlfile | path     | Input mzML file (`*.mzML`)                            |

## Outputs

| Channel       | Type             | Description                            |
| ------------- | ---------------- | -------------------------------------- |
| mzmls_indexed | tuple(val, path) | Indexed mzML file in `out/` directory  |
| versions      | path             | FileConverter version (`versions.yml`) |
| log           | path             | Log file (`*.log`)                     |

## Parameters

| Parameter | Default | Description                                                  |
| --------- | ------- | ------------------------------------------------------------ |
| ext.args  | ''      | Additional CLI arguments appended to `FileConverter` command |

## Container

`community.wave.seqera.io/library/openms:3.6.0--88df48c0a426b7f7` (bioconda OpenMS 3.6.0, the image of nf-core/modules openms modules; `environment.yml`)

## Usage

Used in the `FILE_PREPARATION` subworkflow (`subworkflows/local/file_preparation/main.nf`) for converting mzML files to indexed format before downstream analysis.

## References

- [OpenMS FileConverter](http://www.openms.de/doxygen/nightly/html/TOPP_FileConverter.html)
- [OpenMS homepage](https://www.openms.de/)
