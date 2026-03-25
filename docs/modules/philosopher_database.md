# PHILOSOPHER_DATABASE

Prepares a protein database with decoy sequences and common contaminants using Philosopher.

## Description

This module uses the Philosopher `database` command to add reversed decoy sequences and common contaminant proteins to a FASTA database for proteomics analysis.
It detects whether the input FASTA already contains decoy entries (by searching for the decoy tag prefix) to avoid double-decoys: if decoys exist, it uses `--annotate`; if not, it uses `--custom` to add decoys and contaminants.
Go runtime tuning parameters (GOMEMLIMIT, GOGC, GOMAXPROCS) are set for optimal performance.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| fasta | path | Protein sequence database in FASTA format |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| fasta | tuple(val, path) | `${prefix}_philosopher.fasta` database with decoys/contaminants added |
| log | tuple(val, path) | `philosopher_database.log` command log |
| versions_philosopher | topic: versions | Philosopher version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments for `philosopher database --custom` (only when decoys not present) |
| ext.prefix | `${meta.id}` | Output file prefix |
| ext.decoy_tag | `rev_` | Decoy protein prefix used for detection and generation |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory for Philosopher binary discovery |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_DATABASE` subworkflow (`subworkflows/local/fragpipe_database/main.nf`) as the first step in all FragPipe-based workflows (TMT Label Check, DDA LFQ) to prepare the protein database for downstream search and filtering.

## References

- [Philosopher](https://philosopher.nesvilab.org/)
- [Philosopher wiki](https://github.com/Nesvilab/philosopher/wiki)
- [Publication](https://doi.org/10.1038/s41592-020-0912-y)
