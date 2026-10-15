# PARSEFRAGPIPEWORKFLOW

Parses FragPipe workflow files into per-tool configuration files and unified JSON output.

## Description

PARSEFRAGPIPEWORKFLOW takes a FragPipe `.workflow` file and extracts per-tool configuration into individual config files and a unified JSON output for dynamic workflow orchestration.
It generates configs for all FragPipe tools: msfragger, msbooster, percolator, peptideprophet, proteinprophet, filter, ionquant, tmtintegrator, freequant, diaumpire, diatracer, diann, speclibgen, ptmprophet, ptmshepherd, crystalc, opair, mbg, fpop, saintexpress, skyline, metaproteomics, transferlearning, and database.
The JSON output includes run flags, args, and config_type for each tool, enabling dynamic subworkflow orchestration.

## MSFragger params files

The four MSFragger params files are byte-identical to the ones FragPipe 24.0 writes (`fragger[_dda].params`, `fragger_dia.params`, `fragger_gpfdia.params`, `fragger_dda_plus.params`).
The parser replays FragPipe's own steps: the closed-search defaults (`fragger_closed.params`, copied from FragPipe 24.0) fill the MSFragger tab, the workflow's `msfragger.*` values are loaded in FragPipe's `HashMap` order with the enzyme and deisotope listeners, and the form is written back with `Props.writeProps` (key order, blank lines, disabled mods as `# variable_mod_NN`, `\t\t\t# comment` suffixes).
Keys the workflow does not set get the FragPipe default.
`database_name` is the workflow's `database.db-path` and `num_threads` its `workflow.threads`; MSFRAGGER replaces both values at run time.
`tests/fixtures/fragger_params/` holds the files FragPipe 24.0 wrote for bundled FragPipe workflows (`fragger_params.nf.test`).

## Inputs

| Channel       | Type     | Description                                          |
| ------------- | -------- | ---------------------------------------------------- |
| meta          | val(map) | Sample/experiment metadata map                       |
| workflow_file | path     | FragPipe `.workflow` file containing tool parameters |

## Outputs

| Channel                   | Type                    | Description                                                                        |
| ------------------------- | ----------------------- | ---------------------------------------------------------------------------------- |
| tool_configs_json         | tuple(val, path)        | JSON with all tool configs including run flags and arguments (`tool_configs.json`) |
| msfragger_config          | tuple(val, path)        | MSFragger DDA params file (`configs/fragger.params`)                               |
| msfragger_dia_config      | tuple(val, path)        | MSFragger DIA / DIA-Lib params file (`configs/fragger_dia.params`)                 |
| msfragger_gpfdia_config   | tuple(val, path)        | MSFragger GPF-DIA params file (`configs/fragger_gpfdia.params`)                    |
| msfragger_dda_plus_config | tuple(val, path)        | MSFragger DDA+ params file (`configs/fragger_dda_plus.params`)                     |
| ionquant_config           | tuple(val, path)        | IonQuant config (`configs/ionquant.config`, optional)                              |
| msbooster_config          | tuple(val, path)        | MSBooster params file (`configs/msbooster.config`, optional)                       |
| percolator_config         | tuple(val, path)        | Percolator CLI config (`configs/percolator.config`, optional)                      |
| phi_report_config         | tuple(val, path)        | Philosopher filter/report CLI config (`configs/phi-report.config`, optional)       |
| peptideprophet_config     | tuple(val, path)        | PeptideProphet CLI config (`configs/peptide-prophet.config`, optional)             |
| proteinprophet_config     | tuple(val, path)        | ProteinProphet CLI config (`configs/protein-prophet.config`, optional)             |
| database_config           | tuple(val, path)        | Philosopher database CLI config (`configs/database.config`, optional)              |
| tmtintegrator_config      | tuple(val, path)        | TMTIntegrator YAML config (`configs/tmtintegrator.yml`, optional)                  |
| all_configs               | tuple(val, path)        | Directory containing all generated config files (`configs/`)                       |
| versions_python           | tuple (topic: versions) | Python version                                                                     |

## Parameters

| Parameter | Default | Description                                                                   |
| --------- | ------- | ----------------------------------------------------------------------------- |
| ext.args  | ''      | Passed to the legacy `parse_fragpipe_workflow.py --outdir configs` invocation |

## Container

`community.wave.seqera.io/library/python:3.12.13--27a817c2c0890658`

## Usage

Used in the `FRAGPIPE` subworkflow (`subworkflows/local/fragpipe/main.nf`) to parse `.workflow` files into tool-specific configs that drive downstream module execution.
Also used via the helper functions in `subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils.nf`.

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
