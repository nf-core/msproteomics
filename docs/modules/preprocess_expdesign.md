# PREPROCESS_EXPDESIGN

Preprocesses experimental design files for OpenMS compatibility.

## Description

PREPROCESS_EXPDESIGN converts experimental design file extensions from `.raw` to `.mzML` and extracts a configuration TSV from the design file.
This preprocessing step is needed because OpenMS tools require mzML file references rather than raw vendor formats.
It uses simple `sed`/`grep` commands with the quantms-utils container for environment consistency.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| design | path | Experimental design file in TSV format (`*.tsv` or `*.txt`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| ch_expdesign | path | Design file with `.raw` replaced by `.mzML` (`*_openms_design.tsv`) |
| ch_config | path | Extracted configuration TSV (`*_config.tsv`) |
| versions | path | sdrf-pipelines version (`versions.yml`) |

## Parameters

No configurable parameters. This process uses only `sed`/`grep` commands with no `ext.args` support.

## Container

`biocontainers/quantms-utils:0.0.24--pyh7e72e81_0`

## Usage

Available as a local module for preprocessing experimental design files before OpenMS-based analysis steps.

## References

- [quantms GitHub](https://github.com/bigbio/quantms)
