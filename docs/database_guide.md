# nf-core/msproteomics: Database Preparation

## Overview

Protein identification in mass spectrometry proteomics requires a FASTA database of protein sequences to search spectra against.
The nf-core/msproteomics pipeline supports both pre-built reference proteomes and custom user-provided databases.
Pre-built databases are fetched automatically from UniProt based on organism metadata in the samplesheet, while custom databases can be supplied via the `--database` parameter.

## Pre-Built Reference Proteomes

The pipeline includes pre-built reference proteomes for eight commonly studied organisms.
These databases contain Swiss-Prot reviewed proteins with isoforms, fetched from the [UniProt REST API](https://www.uniprot.org/help/api).

| Organism                   | Common Name | UniProt Proteome ID | Taxon ID |
| -------------------------- | ----------- | ------------------- | -------- |
| _Homo sapiens_             | Human       | UP000005640         | 9606     |
| _Mus musculus_             | Mouse       | UP000000589         | 10090    |
| _Rattus norvegicus_        | Rat         | UP000002494         | 10116    |
| _Saccharomyces cerevisiae_ | Yeast       | UP000002311         | 559292   |
| _Danio rerio_              | Zebrafish   | UP000000437         | 7955     |
| _Drosophila melanogaster_  | Fruit fly   | UP000000803         | 7227     |
| _Caenorhabditis elegans_   | Roundworm   | UP000001940         | 6239     |
| _Escherichia coli_ (K-12)  | E. coli     | UP000000625         | 83333    |

### How pre-built databases are selected

The pipeline auto-selects the correct reference proteome based on the organism metadata in your samplesheet SDRF.
No additional parameters are needed when working with one of the eight supported organisms.
If the organism in your samplesheet matches a pre-built proteome, the pipeline downloads the FASTA automatically at runtime.

You can override the automatic selection by providing `--database /path/to/custom.fasta`.

### Configuration

Pre-built proteome URLs are defined in `conf/reference_proteomes.config`.
Each entry uses the UniProt REST API endpoint:

```
https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&query=(reviewed:true)+AND+(proteome:<PROTEOME_ID>)
```

## Using a Custom Database

To use a custom FASTA database instead of a pre-built reference proteome, provide the path with `--database`:

```bash
nextflow run main.nf --database /path/to/custom.fasta --input samplesheet.csv --outdir results
```

### Supported formats

The pipeline accepts standard FASTA files with any of the following header formats:

- **UniProt**: `>sp|P12345|PROTEIN_NAME OS=Homo sapiens`
- **NCBI**: `>gi|123456|ref|NP_000001.1| protein description`
- **Generic**: `>protein_id description`

### Downloading from UniProt

You can download a Swiss-Prot reviewed database with isoforms directly from the UniProt REST API.

Example for human:

```bash
curl -o human_swissprot.fasta \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&query=(reviewed:true)+AND+(proteome:UP000005640)"
```

Alternatively, use the [UniProt web interface](https://www.uniprot.org/) to search and download sequences with custom filters.

## Contaminant Sequences

Contaminant proteins (keratins, trypsin, BSA, and other common lab contaminants) are important to include in the search database to prevent false identifications.

### FragPipe mode (`--mode fragpipe`)

Contaminants are added **automatically** by Philosopher's database command.
No user action is required.
Philosopher appends a standard set of common lab contaminant sequences to the database during the database preparation step.

### DIA-NN mode (`--mode diann`)

Contaminants are **NOT** added automatically.
Users should append contaminant sequences to their FASTA file before running the pipeline.

The [cRAP (common Repository of Adventitious Proteins)](https://www.thegpm.org/crap/) database is the standard contaminant collection for proteomics.

To append cRAP contaminants to your database:

```bash
curl -o crap.fasta https://www.thegpm.org/crap/crap.fasta
cat your_database.fasta crap.fasta > database_with_contaminants.fasta
```

Then provide the combined file:

```bash
nextflow run main.nf --mode diann --database database_with_contaminants.fasta --input samplesheet.csv --outdir results
```

## Decoy Sequences

Decoy sequences (typically reversed protein sequences) are required for target-decoy FDR estimation.
Users do **not** need to add decoys manually for either pipeline mode.

- **FragPipe mode**: Decoy sequences (reversed, with `rev_` prefix) are generated automatically by Philosopher's database command.
- **DIA-NN mode**: Decoy sequences are generated internally by DIA-NN during analysis.

## Best Practices

- **Swiss-Prot vs TrEMBL**: Swiss-Prot (reviewed) sequences are recommended for most analyses because they are manually curated and non-redundant.
  TrEMBL (unreviewed) sequences may be needed for non-model organisms with limited Swiss-Prot coverage.
- **Isoforms**: Include isoforms (`includeIsoform=true`) for comprehensive protein coverage, especially for quantitative studies where isoform-level resolution matters.
- **Database size**: Smaller, well-curated databases improve statistical sensitivity by reducing the multiple testing burden.
  Avoid using the entire UniProt KB unless your experiment specifically requires it.
- **Species matching**: Ensure the database species matches your sample species.
  Using a mismatched database will result in poor identification rates.
- **Version tracking**: Record the UniProt release date used for each analysis to ensure reproducibility.
  UniProt releases new versions approximately every 8 weeks.
- **Non-model organisms**: For organisms not in the pre-built list, download a species-specific FASTA from UniProt and provide it via `--database`.

## Organism Name Resolution

The pipeline recognizes common organism name formats and normalizes them automatically.
This mapping is handled by `WorkflowUtils.convertOrganismToStandardName()`.

The following aliases are currently supported:

| Standard Name              | Accepted Aliases                          |
| -------------------------- | ----------------------------------------- |
| _Homo sapiens_             | `human`, `hs`, `Homo_sapiens`             |
| _Mus musculus_             | `mouse`, `mm`, `Mus_musculus`             |
| _Saccharomyces cerevisiae_ | `yeast`, `sc`, `Saccharomyces_cerevisiae` |

For other supported organisms, use the full binomial name (e.g., `Rattus norvegicus`, `Danio rerio`) in your samplesheet.
Underscores and spaces are treated equivalently, and matching is case-insensitive.
