# MERGESPLITSEARCH

Merge MSFragger split-database search results per sample.

## Description

MERGESPLITSEARCH merges pepXML, PIN, and score histogram files from all database chunks for a single sample, sums histograms, generates expect functions via MSFragger, and re-ranks hits per spectrum.
The merge logic is faithfully ported from FragPipe's `msfragger_pep_split.py`, sorting search hits by (expectscore, 1/hyperscore, abs(massdiff)) and keeping the top N results.
Used within the FRAGPIPE_SEARCH subworkflow when MSFragger runs with split-database mode.

## Inputs

| Channel          | Type     | Description                                                                                                                                            |
| ---------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| meta             | val(map) | Sample metadata map (e.g., `[ id:'sample1' ]`)                                                                                                         |
| chunk_pepxmls    | path     | pepXML files from all database chunks for this sample (chunk-prefixed)                                                                                 |
| chunk_pins       | path     | PIN (Percolator input) files from all database chunks (chunk-prefixed)                                                                                 |
| chunk_histograms | path     | Score histogram TSV files from all database chunks (chunk-prefixed)                                                                                    |
| params_file      | path     | The file's search params; `output_report_topN`, `output_max_expect` and `output_file_extension` are read from it (`msfragger_pep_split.py:76-80`)      |
| num_chunks       | val(int) | Number of database chunks used in the split search                                                                                                     |
| fasta            | path     | The full FASTA; its path in the task (`${PWD}/<fasta>`) replaces the chunk FASTA in the merged pepXML header (`msfragger_pep_split.py:92-94, 302-303`) |
| msfragger_dir    | path     | Unzipped MSFragger tool directory (optional, pass `[]` when not using)                                                                                 |

## Outputs

| Channel            | Type              | Description                                                                    |
| ------------------ | ----------------- | ------------------------------------------------------------------------------ |
| pepxml             | tuple(meta, path) | Merged pepXML file with re-ranked search hits (`<prefix>/*.pepXML`)            |
| pin                | tuple(meta, path) | Merged PIN file with re-ranked and reconciled protein lists (`<prefix>/*.pin`) |
| license_agreement  | path              | License agreement marker file                                                  |
| versions_python    | tuple             | Python software version (topic channel)                                        |
| versions_msfragger | tuple             | MSFragger software version (topic channel)                                     |
| versions_fragpipe  | tuple             | FragPipe software version (topic channel)                                      |

## Parameters

| Parameter                            | Default         | Description                                       |
| ------------------------------------ | --------------- | ------------------------------------------------- |
| ext.prefix                           | meta.id         | Override for the sample name                      |
| ext.fragpipe_tools_dir               | (auto-detected) | Override for the FragPipe tools directory path    |
| ext.java_xmx                         | (auto)          | Override for Java heap size                       |
| ext.agree_fragpipe_license_agreement | (required)      | Must be `true` or the process exits with an error |

## Container

The licensed FragPipe image given with `--fragpipe_container`.

## Usage

Used in the `FRAGPIPE_SEARCH` subworkflow (`subworkflows/local/fragpipe_search/main.nf`) to merge results when MSFragger uses split-database search mode for large FASTA databases.

## References

- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)
- [FragPipe homepage](https://fragpipe.nesvilab.org/)
