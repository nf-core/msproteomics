# PMULTIQC

Proteomics-specific QC report generation using the pmultiqc MultiQC plugin.

## Description

PMULTIQC generates proteomics quality control reports using the pmultiqc library, which extends the MultiQC framework with proteomics-specific modules.
It produces an HTML report with optional database output, supporting configurable quantification methods, contaminant filtering, and table/idxml plot toggling.
This module aggregates results from upstream pipeline stages into a comprehensive QC report.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| results | path | All result files staged into `results/` directory |
| quantms_log | path | quantms log file |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| ch_pmultiqc_report | path | MultiQC HTML report (`*.html`) |
| ch_pmultiqc_db | path | SQLite3 database with protein, PSM, and quantification data (`*.db`, optional) |
| versions | path | pmultiqc version (`versions.yml`) |
| data | path | MultiQC data directory (`*_data`) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments passed to the `multiqc` command |
| params.enable_pmultiqc | - | Combined with `params.export_mztab` to control `--quantms_plugin` |
| params.skip_table_plots | - | Controls `--disable_table` flag |
| params.pmultiqc_idxml_skip | - | Controls `--ignored_idxml` flag |
| params.contaminant_string | - | Sets `--contaminant_affix` for contaminant filtering |
| params.quantification_method | - | Sets `--quantification_method` for report generation |

## Container

`biocontainers/pmultiqc:0.0.39--pyhdfd78af_0`

## Usage

Available as a local module for generating proteomics QC reports.
Designed to aggregate outputs from DIA and DDA workflows into a single HTML report.

## References

- [pmultiqc GitHub](https://github.com/bigbio/pmultiqc/)
- [MultiQC homepage](https://multiqc.info/)
