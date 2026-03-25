# METAPROTEOMICS

Metaproteomics database optimization and taxonomy analysis.

## Description

The Metaproteomics module performs database optimization for metaproteomics workflows using FragPipe's FP-Meta DbOptimizer tool.
It analyzes FragPipe search results to identify taxonomic composition and generates an optimized FASTA database by filtering proteins based on taxonomic assignments, enabling iterative refinement of the search database.
Used within the FRAGPIPE_EXPORT subworkflow for metaproteomics experiments requiring NCBI taxonomy files.

## Inputs

| Channel | Type | Description |
|---------|------|-------------|
| meta | val(map) | Sample metadata map |
| project_dir | path | Directory containing FragPipe results for taxonomy analysis |
| config_cli | val(string) | Metaproteomics parameters as CLI string (e.g., `--decoyTag rev_ --qvalue 0.01 --iterations 2`) |
| fasta | path | Protein sequence database in FASTA format (shared resource) |
| taxon_name_file | path | NCBI taxonomy `names.dmp` file (shared resource) |
| taxon_node_file | path | NCBI taxonomy `nodes.dmp` file (shared resource) |

## Outputs

| Channel | Type | Description |
|---------|------|-------------|
| results_dir | tuple(meta, path) | Directory containing metaproteomics results and log |
| optimized_fasta | tuple(meta, path) | Optimized FASTA database filtered by taxonomic analysis (`*_optimized.fasta`, optional) |
| taxonomy_results | tuple(meta, path) | Taxonomy analysis results (`*_taxonomy*.tsv`, optional) |
| results | tuple(meta, path) | Metaproteomics analysis results (`*_metaproteomics*.tsv`, optional) |
| versions_metaproteomics | tuple | Metaproteomics (FP-Meta) software version (topic channel) |
| versions_fragpipe | tuple | FragPipe software version (topic channel) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ext.args | '' | Additional CLI arguments appended to the FP-Meta DbOptimizer Java command |
| ext.prefix | meta.id | Override for the output directory name |
| ext.fragpipe_tools_dir | (auto-detected) | Override for the FragPipe tools directory path |

## Container

`docker.io/fcyucn/fragpipe:24.0`

## Usage

Used in the `FRAGPIPE_EXPORT` subworkflow (`subworkflows/local/fragpipe_export/main.nf`) for metaproteomics database optimization.
Enabled conditionally via `shouldRunTool(configs, 'metaproteomics')`.
Requires NCBI taxonomy files (`names.dmp` and `nodes.dmp`) as shared resources.

## References

- [FragPipe homepage](https://fragpipe.nesvilab.org/)
- [Metaproteomics tutorial](https://fragpipe.nesvilab.org/docs/tutorial_metaproteomics.html)
- [FragPipe GitHub](https://github.com/Nesvilab/FragPipe)
