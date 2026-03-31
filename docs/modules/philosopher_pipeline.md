# PHILOSOPHER_PIPELINE

Combined Philosopher filter and report pipeline in a single process.

## Description

This module runs the complete Philosopher filter and report workflow including workspace initialization, database annotation, FDR filtering, and result reporting in a single process.
It is similar to PHILOSOPHER_FILTER but also accepts mzML files as input and always runs the report step.
All generated TSV files are moved to the output directory.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| pepxml_files | path | PepXML files from search/validation (`*.pep.xml`) |
| protxml | path | ProtXML from ProteinProphet (`*.prot.xml`) |
| mzml_files | path | mzML spectral files |
| filter_flags | val(string) | Filter CLI flags (e.g., `--sequential --razor --prot 0.01 --picked`) |
| report_flags | val(string) | Report CLI flags (e.g., `--msstats`) |
| fasta | path | Philosopher database with decoys/contaminants from PHILOSOPHER_DATABASE |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Output directory `${prefix}` containing all result files |
| proteins | tuple(val, path) | `${prefix}/protein.tsv` protein-level results |
| peptides | tuple(val, path) | `${prefix}/peptide.tsv` peptide-level results |
| psms | tuple(val, path) | `${prefix}/psm.tsv` PSM-level results |
| versions_philosopher | topic: versions | Philosopher version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments for `philosopher filter` |
| ext.prefix | `${meta.id}` | Output directory name |
| ext.decoy_tag | `rev_` | Decoy protein prefix |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Available as a convenience process that combines PHILOSOPHER_FILTER with report into a single step, useful when mzML files are needed alongside pepXML/protXML for processing.

## References

- [Philosopher](https://philosopher.nesvilab.org/)
- [Philosopher wiki](https://github.com/Nesvilab/philosopher/wiki)
- [Publication](https://doi.org/10.1038/s41592-020-0912-y)
