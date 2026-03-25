# FRAGPIPE

Run complete FragPipe proteomics workflow as an all-in-one process.

## Description

FragPipe is a comprehensive proteomics analysis suite that orchestrates MSFragger, Philosopher, IonQuant, and other tools for DDA/DIA/TMT analysis.
This module runs the entire FragPipe pipeline using the headless CLI, copying input files locally for mzBIN cache writes, creating a manifest file, and launching with full tool configuration.
Used in the main `msproteomics.nf` workflow and the FRAGPIPE_WF subworkflow for end-to-end proteomics analysis.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map (e.g., `[ id:'experiment1' ]`) |
| mzml_files | path | Mass spectrometry data files in mzML format |
| fasta | path | Protein sequence database in FASTA format |
| workflow_file | path | FragPipe `.workflow` configuration file |
| data_type | val(string) | Data type: `DDA`, `DIA`, `GPF-DIA`, `DIA-Quant`, `DIA-Lib`, `DDA+` |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(meta, path) | Complete FragPipe output directory |
| proteins | tuple(meta, path) | Protein-level quantification results (`*/protein.tsv`) |
| peptides | tuple(meta, path) | Peptide-level quantification results (`*/peptide.tsv`) |
| psms | tuple(meta, path) | Peptide-spectrum match results (`*/psm.tsv`) |
| ions | tuple(meta, path) | Ion reports (`*/ion.tsv`, optional) |
| pepxml | tuple(meta, path) | pepXML search results (`*/*.pepXML`, optional) |
| protxml | tuple(meta, path) | Combined ProtXML (`combined.prot.xml`, optional) |
| versions_fragpipe | tuple | FragPipe software version (topic channel) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the `fragpipe --headless` command |
| ext.prefix | meta.id | Override for the output directory name |

## Container

Inherited from pipeline configuration (requires the full FragPipe container).

## Usage

Used in the main workflow (`workflows/msproteomics.nf`) as `FRAGPIPE_PIPELINE` and `FRAGPIPE_PIPELINE_TMT` via the `FRAGPIPE_WF` subworkflow.
Supports all FragPipe workflow types including DDA, DIA, TMT, and label-free quantification.

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
- [FragPipe tutorial](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe.html)
- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)
- [Publication: doi:10.1038/s41592-020-0967-5](https://doi.org/10.1038/s41592-020-0967-5)
