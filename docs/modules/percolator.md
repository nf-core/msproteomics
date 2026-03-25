# PERCOLATOR

Semi-supervised PSM validation and FDR control using Percolator.

## Description

Percolator is a semi-supervised machine learning tool that uses support vector machines to discriminate between correct and incorrect peptide-spectrum matches (PSMs), enabling accurate FDR estimation.
This module runs Percolator followed by conversion to pepXML format using FragPipe's `PercolatorOutputToPepXML` converter.
It supports both PSM-level (`--only-psms`) and peptide-level output modes and auto-detects DDA+ ion mobility data (ranked pepXMLs).

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| pin_file | path | Percolator input (PIN) file from MSFragger/MSBooster (`*.pin`) |
| pepxml_files | path | PepXML files from MSFragger (for perc2pepxml conversion, `*.pepXML`) |
| mzml_file | path | Mass spectrometry data file (for perc2pepxml path, `*.mzML`) |
| percolator_cli | val(string) | CLI flags for Percolator (e.g., `--only-psms --no-terminate --post-processing-tdc`) |
| data_type | val(string) | `'DDA'` or `'DIA'` - affects pepXML output format |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(val, path) | Results directory containing all Percolator output files and logs |
| pepxml | tuple(val, path) | Converted pepXML results with probabilities (`interact-*.pep.xml`) |
| target_psms | tuple(val, path) | Target PSM results (`*_target_psms.tsv`) |
| decoy_psms | tuple(val, path) | Decoy PSM results (`*_decoy_psms.tsv`) |
| versions_percolator | tuple (topic: versions) | Percolator version |
| versions_fragpipe | tuple (topic: versions) | FragPipe version |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to both Percolator and pepXML conversion commands |
| ext.min_prob | 0.5 | Minimum probability for pepXML conversion |
| ext.decoy_tag | (default) | Sets `--protein-decoy-pattern` for Percolator |
| ext.fragpipe_tools_dir | (default path) | Override FragPipe tools directory |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_VALIDATE` subworkflow (`subworkflows/local/fragpipe_validate/main.nf`) for PSM-level FDR control.
The module auto-detects `--only-psms` in the combined `percolator_cli + args` to choose PSM-level vs peptide-level output flags.
Auto-switches to DIA mode for the pepXML converter when ranked pepXMLs are detected (DDA+ ion mobility).

## References

- [Percolator homepage](http://percolator.ms/)
- [Percolator Wiki](https://github.com/percolator/percolator/wiki)
- [Publication: Kall et al., Nature Methods, 2007](https://doi.org/10.1038/nmeth1113)
