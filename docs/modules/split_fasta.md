# SPLIT_FASTA

Splits a FASTA database into M chunks at protein boundaries for parallel MSFragger search.

## Description

This module splits a FASTA file into M chunks at protein boundaries using memory-mapped I/O via the `split_fasta.py` bin script.
It replicates `get_fasta_offsets()` and `set_up_directories()` from FragPipe's `msfragger_pep_split.py` for split-database search.
Each chunk is written to a numbered subdirectory under `split_db/`.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| fasta | path | Protein sequence database in FASTA format |
| num_chunks | val(integer) | Number of chunks to split into (>= 1) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| fasta_chunks | tuple(val, path) | `split_db/*/*.fasta` chunk files in numbered subdirectories |
| versions_python | topic: versions | Python version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional arguments (not directly used by the bin script) |
| ext.prefix | `${meta.id}` | Prefix for naming |

## Container

`quay.io/biocontainers/python:3.12`

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`) to split the FASTA database for parallel MSFragger split-database search, improving performance on large databases.

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
- [Split database search tutorial](https://fragpipe.nesvilab.org/docs/tutorial_split_db.html)
