# PTMPROPHET

PTM site localization using PTMProphet from the Trans-Proteomic Pipeline.

## Description

PTMProphet computes probabilities for the precise localization of post-translational modifications on peptide sequences using pepXML search results as input.
It is a tool from the Trans-Proteomic Pipeline (TPP), bundled in the FragPipe distribution.
The module always uses `MAXTHREADS=1` (matching FragPipe convention, since Nextflow handles parallelism at the process level).

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| pepxml_file | path | PepXML file containing peptide search results with PTMs (`*.pep.xml`) |
| config_cli | val(string) | PTMProphet parameters as CLI string (UPPERCASE=value format, e.g., `MINPROB=0.5 STATIC KEEPOLD`) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Results directory containing PTMProphet output and logs |
| mod_pepxml | tuple(val, path) | Modified pepXML with PTM site localization probabilities (`*.mod.pep.xml`) |
| versions_ptmprophet | tuple (topic: versions) | PTMProphet version |
| versions_fragpipe | tuple (topic: versions) | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Inserted into PTMProphet command after `config_cli` and before input/output files |
| ext.decoy_tag | 'rev_' | Decoy prefix for decoy protein identification |
| ext.fragpipe_tools_dir | (default path) | Override FragPipe tools directory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_VALIDATE` subworkflow (`subworkflows/local/fragpipe_validate/main.nf`) for PTM site localization after Percolator validation.
Any `MAXTHREADS=` value in `config_cli` is stripped and replaced with `MAXTHREADS=1`.
The module copies the input pepXML to the output directory, runs PTMProphet, then cleans up the copy.

## References

- [Trans-Proteomic Pipeline (TPP)](http://tools.proteomecenter.org/wiki/index.php?title=Software:TPP)
- [TPP SourceForge](https://sourceforge.net/projects/sashimi/)
