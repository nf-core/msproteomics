# SAMPLESHEET_CHECK

Validates experimental design files (SDRF or samplesheet) using quantms-utils.

## Description

This module validates samplesheet or SDRF-format experimental design files using `quantmsutilsc checksamplesheet`.
It supports optional skipping of SDRF validation, MS validation, factor validation, experimental design validation, and OLS cache-only mode via pipeline parameters.
The validated file is passed through to downstream processes.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| input_file | path | Samplesheet or SDRF file (`*.tsv`, `*.txt`, `*.csv`) |
| is_sdrf | val(boolean) | Whether the input file is in SDRF format |
| validate_ontologies | val(boolean) | Whether to validate ontologies (false skips SDRF validation) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| log | path | `*.log` validation log |
| checked_file | path | Pass-through of the validated input file |
| versions | path | `versions.yml` with software versions |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments for `quantmsutilsc checksamplesheet` |
| params.skip_ms_validation | false | Skip mass spectrometry validation |
| params.skip_factor_validation | false | Skip factor validation |
| params.skip_experimental_design_validation | false | Skip experimental design validation |
| params.use_ols_cache_only | false | Use OLS cache only for ontology lookups |

## Container

`biocontainers/quantms-utils:0.0.24--pyh7e72e81_0`

## Usage

Used in the `INPUT_CHECK` subworkflow (`subworkflows/local/input_check/main.nf`) as the first step to validate user-provided experimental design files before pipeline execution.

## References

- [quantms-utils](https://github.com/bigbio/quantms)
