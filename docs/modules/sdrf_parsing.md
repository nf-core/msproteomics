# SDRF_PARSING

Converts SDRF proteomics files into OpenMS experimental design format.

## Description

This module parses an SDRF file using `parse_sdrf convert-openms` from the sdrf-pipelines package.
It produces a configuration TSV and an experimental design TSV suitable for downstream analysis tools.
When `params.convert_dotd` is enabled, it handles `.d` format extension conversions to mzML.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| sdrf | path | A valid SDRF file |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| ch_expdesign | path | `${sdrf.baseName}_openms_design.tsv` experimental design file in OpenMS format |
| ch_sdrf_config_file | path | `${sdrf.baseName}_config.tsv` config file with search engine parameters |
| log | path | `*.log` parsing log |
| versions | path | `versions.yml` with software versions |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments appended to `parse_sdrf convert-openms` |
| params.convert_dotd | false | Enable `.d` format extension conversions to mzML |

## Container

`biocontainers/sdrf-pipelines:0.0.33--pyhdfd78af_0`

## Usage

Used in the `CREATE_INPUT_CHANNEL` subworkflow (`subworkflows/local/create_input_channel/main.nf`) to parse SDRF files into configuration and experimental design TSVs for downstream processing.

## References

- [sdrf-pipelines](https://github.com/bigbio/sdrf-pipelines)
