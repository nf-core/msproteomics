# GENERATE_SDRF_FROM_SAMPLESHEET

Convert a samplesheet CSV to SDRF (Sample and Data Relationship Format) TSV.

## Description

GENERATE_SDRF_FROM_SAMPLESHEET converts a simple CSV samplesheet into an SDRF-proteomics TSV file using the `generate_sdrf.py` script.
It takes pipeline parameters as a JSON string to populate SDRF columns including organism, enzyme, modifications, and other experiment-level metadata.
Used at the entry point of the main workflow when the user provides a samplesheet instead of a pre-built SDRF file.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| samplesheet | path | Input samplesheet CSV with columns: sample, spectra, condition (and optional label, fraction, replicate) |
| params_json | val(string) | Pipeline parameters as JSON string containing experiment-level parameters |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| sdrf | path | Generated SDRF-proteomics file (`*.sdrf.tsv`) |
| versions | path | Software versions (`versions.yml`) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the `generate_sdrf.py` command |

## Container

`biocontainers/python:3.12`

## Usage

Used in the main workflow (`workflows/msproteomics.nf`) to generate SDRF files from user-provided samplesheets.
This enables users to provide a simple CSV format instead of the full SDRF specification.

## References

- [SDRF-Proteomics specification](https://github.com/bigbio/proteomics-sample-metadata)
- [Python](https://www.python.org)
