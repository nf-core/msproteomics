# PROTEINPROPHET

Statistical protein inference using ProteinProphet via Philosopher.

## Description

ProteinProphet is a statistical model for protein inference from peptide identifications, handling shared peptides and computing protein-level probabilities considering peptide degeneracy.
This module takes multiple pepXML files, writes them to a filelist, initializes a Philosopher workspace, and runs ProteinProphet with `--output combined`.
It produces `combined.prot.xml` and passes through the input pepXML files for downstream use.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample/experiment metadata |
| pepxml_files | path | PepXML files from Percolator or PeptideProphet (`*.pep.xml`) |
| config_cli | val(string) | ProteinProphet CLI options (e.g., `--maxppmdiff 2000000`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Output directory `${prefix}` with ProteinProphet results and logs |
| protxml | tuple(val, path) | `${prefix}/combined.prot.xml` protein inference results |
| pepxml | tuple(val, path) | Pass-through of input pepXML files for downstream synchronization |
| versions_philosopher | topic: versions | Philosopher version |
| versions_fragpipe | topic: versions | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional command-line arguments appended after `config_cli` and `--output combined` |
| ext.prefix | `${meta.id}` | Output directory name |
| ext.fragpipe_tools_dir | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools` | FragPipe tools directory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_INFERENCE` subworkflow (`subworkflows/local/fragpipe_inference/main.nf`) for aggregate protein inference across all samples before per-sample PHILOSOPHER_FILTER.

## References

- [ProteinProphet](https://proteinprophet.sourceforge.net/)
- [Publication](https://doi.org/10.1021/ac0341261)
