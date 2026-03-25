# TMTINTEGRATOR

Runs TMT-Integrator for isobaric labeling quantification.

## Description

TMT-Integrator extracts and combines channel abundances from multiple TMT-labeled samples, producing comprehensive quantification reports at gene, protein, peptide, and PTM site levels.
It supports both single-plex and multi-plex modes, with automatic detection based on annotation file format (6+ column TMTIntegrator format vs 2-column IonQuant format).
The module generates abundance and ratio tables and handles annotation format conversion between IonQuant and TMTIntegrator formats.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| psm_dirs | path (stageAs `input_results/*`) | Per-sample directories containing psm.tsv and protein.tsv from PHILOSOPHER_FILTER |
| annotation_files | path | TMT annotation file (2-column IonQuant format or 6-column TMTIntegrator format) |
| config_file | path | TMT-Integrator YAML configuration file (ref_tag, min_pep_prob, prot_norm, groupby, etc.) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Output directory `${prefix}` containing all TMT-Integrator results |
| abundance | tuple(val, path) | `${prefix}/abundance_*.tsv` abundance tables at various levels (optional) |
| ratio | tuple(val, path) | `${prefix}/ratio_*.tsv` ratio tables with log2 fold changes (optional) |
| proteins_txt | tuple(val, path) | `${prefix}/*_Proteins.txt` protein-level summary (optional) |
| peptides_txt | tuple(val, path) | `${prefix}/*_Peptides.txt` peptide-level summary (optional) |
| versions_tmtintegrator | topic: versions | TMT-Integrator version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments appended to `java -jar TMT-Integrator.jar` |
| ext.prefix | `${meta.id}` | Output directory name |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory for JAR discovery |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_QUANT` subworkflow (`subworkflows/local/fragpipe_quant/main.nf`) for TMT isobaric labeling quantification after PHILOSOPHER_LABELQUANT and PHILOSOPHER_REPORT have generated per-sample results with reporter ion intensities.

## References

- [TMT-Integrator](https://github.com/Nesvilab/TMT-Integrator)
- [TMT tutorial](https://fragpipe.nesvilab.org/docs/tutorial_tmt.html)
