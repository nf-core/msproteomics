# MERGE_SPLIT_SEARCH

Merge MSFragger split-database search results per sample.

## Description

MERGE_SPLIT_SEARCH merges pepXML, PIN, and score histogram files from all database chunks for a single sample, sums histograms, generates expect functions via MSFragger, and re-ranks hits per spectrum.
The merge logic is faithfully ported from FragPipe's `msfragger_pep_split.py`, sorting search hits by (expectscore, 1/hyperscore, abs(massdiff)) and keeping the top N results.
Used within the FRAGPIPE_SEARCH subworkflow when MSFragger runs with split-database mode.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map (e.g., `[ id:'sample1' ]`) |
| chunk_pepxmls | path | pepXML files from all database chunks for this sample (chunk-prefixed) |
| chunk_pins | path | PIN (Percolator input) files from all database chunks (chunk-prefixed) |
| chunk_histograms | path | Score histogram TSV files from all database chunks (chunk-prefixed) |
| num_chunks | val(int) | Number of database chunks used in the split search |
| params_file | path | MSFragger params file (for `output_report_topN`, `output_max_expect` settings) |
| msfragger_dir | path | Unzipped MSFragger tool directory (optional, pass `[]` when not using) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| pepxml | tuple(meta, path) | Merged pepXML file with re-ranked search hits (`merged/*.pepXML`) |
| pin | tuple(meta, path) | Merged PIN file with re-ranked and reconciled protein lists (`merged/*.pin`) |
| license_agreement | path | License agreement marker file |
| versions_python | tuple | Python software version (topic channel) |
| versions_msfragger | tuple | MSFragger software version (topic channel) |
| versions_fragpipe | tuple | FragPipe software version (topic channel) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.prefix | meta.id | Override for the sample name |
| ext.output_report_topN | 1 | Number of top hits to keep per spectrum |
| ext.output_max_expect | 50.0 | Maximum expect value threshold for merged results |
| ext.fasta_path | (none) | FASTA path for MSFragger expect function generation |
| ext.fragpipe_tools_dir | (auto-detected) | Override for the FragPipe tools directory path |
| ext.java_xmx | (auto) | Override for Java heap size |
| ext.agree_fragpipe_license_agreement | (required) | Must be `true` or the process exits with an error |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`) to merge results when MSFragger uses split-database search mode for large FASTA databases.

## References

- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)
- [FragPipe homepage](https://fragpipe.nesvilab.org/)
