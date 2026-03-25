# DECOMPRESS

Extract compressed archive files into Bruker `.d` directories.

## Description

DECOMPRESS extracts compressed archive files (tar.gz, gz, tar, zip) into `.d` directories (Bruker timsTOF format).
It includes tar verification after extraction and uses a retry strategy with escalating stageInMode (symlink to copy on AWS Batch, link to symlink to copy locally).
Used within the FILE_PREPARATION subworkflow for handling compressed Bruker raw files.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map (uses `meta.mzml_id`) |
| compressed_file | path | Compressed archive file (`.tar.gz`, `.gz`, `.tar`, `.zip`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| decompressed_files | tuple(meta, path) | Decompressed `.d` directory |
| versions | path | Software versions (`versions.yml`) for gunzip, tar, unzip |
| log | path | Decompression log file |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.prefix | meta.mzml_id | Override for the log file prefix |

## Container

`quay.io/biocontainers/mulled-v2-796b0610595ad1995b121d0b85375902097b78d4:a3a3220eb9ee55710d743438b2ab9092867c98c6-0`

## Usage

Used in the `FILE_PREPARATION` subworkflow (`subworkflows/local/file_preparation/main.nf`) to decompress Bruker `.d` archives before downstream processing.
Compressed files are branched and routed through DECOMPRESS, then mixed back with uncompressed files.

## References

- [GNU tar](https://www.gnu.org/software/tar/)
- [GNU gzip](https://www.gnu.org/software/gzip/)
