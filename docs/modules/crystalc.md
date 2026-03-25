# CRYSTALC

Remove chimeric artifact PSMs from pepXML search results.

## Description

Crystal-C is a post-processing tool for tandem mass spectrometry database search results that identifies and removes chimeric spectrum artifacts.
It improves the accuracy of peptide identifications by filtering out PSMs that originate from co-fragmented precursor ions.
Used within the FRAGPIPE_SEARCH subworkflow after MSFragger database search.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| pepxml_file | path | Search results in pepXML format from MSFragger (`*.pepXML`) |
| mzml_file | path | Mass spectrometry data file for spectrum access (`*.mzML`, `*.mzXML`) |
| config_file | path | Crystal-C parameters file with paths configured (`*.params`) |
| fasta | path | Protein sequence database FASTA (shared resource) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(meta, path) | Directory containing Crystal-C results and log |
| pepxml_filtered | tuple(meta, path) | Filtered pepXML files with chimeric artifacts removed (`*_c.pepXML`) |
| versions_crystalc | tuple | Crystal-C software version (topic channel) |
| versions_fragpipe | tuple | FragPipe software version (topic channel) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the `crystalc.Run` Java command |
| ext.prefix | meta.id | Override for the output directory name |
| ext.fragpipe_tools_dir | (auto-detected) | Override for the FragPipe tools directory path |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`) to filter chimeric artifacts from MSFragger search results.
Crystal-C runs per-sample after MSFragger and before validation (Percolator/PeptideProphet).
Enabled conditionally via `shouldRunTool(configs, 'crystalc')`.

## References

- [Crystal-C GitHub](https://github.com/AimeeD90/Crystal-C)
