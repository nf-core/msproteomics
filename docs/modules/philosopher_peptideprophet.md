# PEPTIDEPROPHET

Statistical validation of peptide assignments using PeptideProphet via Philosopher.

## Description

PeptideProphet is a statistical model for peptide-spectrum match (PSM) validation.
It computes probabilities for peptide assignments based on search engine scores and assigns confidence levels to identifications.
This module parses a config file for `peptideprophet=` flags, initializes a Philosopher workspace, runs PeptideProphet with the database and decoy tag, then normalizes `base_name` attributes in the output pepXML to remove absolute Nextflow work directory paths.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| pepxml | path | PepXML file from search engine (`*.pepXML`) |
| mzml | path | mzML spectral file (used for base_name normalization) |
| config_file | path | Config file with `peptideprophet=flags` line from PARSE_FRAGPIPE_WORKFLOW |
| fasta | path | Protein sequence database |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Output directory `${prefix}` with PeptideProphet results and logs |
| pepxml | tuple(val, path) | `${prefix}/interact-*.pep.xml` PepXML with probabilities |
| versions_philosopher | topic: versions | Philosopher version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments for `philosopher peptideprophet` |
| ext.prefix | `${meta.id}` | Output directory name |
| ext.decoy_tag | `rev_` | Decoy protein prefix |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_VALIDATE` subworkflow (`subworkflows/local/fragpipe_validate/main.nf`) as an optional alternative to Percolator for statistical PSM validation.

## References

- [PeptideProphet](https://peptideprophet.sourceforge.net/)
- [Publication](https://doi.org/10.1021/ac025747h)
