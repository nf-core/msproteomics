# PARSE_FRAGPIPE_WORKFLOW

Parses FragPipe workflow files into per-tool configuration files and unified JSON output.

## Description

PARSE_FRAGPIPE_WORKFLOW takes a FragPipe `.workflow` file and extracts per-tool configuration into individual config files and a unified JSON output for dynamic workflow orchestration.
It generates configs for all FragPipe tools: msfragger, msbooster, percolator, peptideprophet, proteinprophet, filter, ionquant, tmtintegrator, freequant, diaumpire, diatracer, diann, speclibgen, ptmprophet, ptmshepherd, crystalc, opair, mbg, fpop, saintexpress, skyline, metaproteomics, transferlearning, and database.
The JSON output includes run flags, args, and config_type for each tool, enabling dynamic subworkflow orchestration.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata map |
| workflow_file | path | FragPipe `.workflow` file containing tool parameters |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| tool_configs_json | tuple(val, path) | JSON with all tool configs including run flags and arguments (`tool_configs.json`) |
| msfragger_config | tuple(val, path) | MSFragger params file (`configs/fragger.params`, optional) |
| ionquant_config | tuple(val, path) | IonQuant config (`configs/ionquant.config`, optional) |
| msbooster_config | tuple(val, path) | MSBooster params file (`configs/msbooster.config`, optional) |
| percolator_config | tuple(val, path) | Percolator CLI config (`configs/percolator.config`, optional) |
| phi_report_config | tuple(val, path) | Philosopher filter/report CLI config (`configs/phi-report.config`, optional) |
| peptideprophet_config | tuple(val, path) | PeptideProphet CLI config (`configs/peptide-prophet.config`, optional) |
| proteinprophet_config | tuple(val, path) | ProteinProphet CLI config (`configs/protein-prophet.config`, optional) |
| database_config | tuple(val, path) | Philosopher database CLI config (`configs/database.config`, optional) |
| tmtintegrator_config | tuple(val, path) | TMTIntegrator YAML config (`configs/tmtintegrator.yml`, optional) |
| all_configs | tuple(val, path) | Directory containing all generated config files (`configs/`) |
| versions_python | tuple (topic: versions) | Python version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Passed to the legacy `parse_fragpipe_workflow.py --outdir configs` invocation |

## Container

`quay.io/biocontainers/python:3.12`

## Usage

Used in the `FRAGPIPE` subworkflow (`subworkflows/local/fragpipe/main.nf`) and `TMT_LABELCHECK` subworkflow (`subworkflows/local/tmt_labelcheck/main.nf`) to parse `.workflow` files into tool-specific configs that drive downstream module execution.
Also used via the `fragpipe_utils.nf` helper subworkflow.

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
- [nf-fragpipe GitHub](https://github.com/nf-fragpipe/nf-fragpipe)
