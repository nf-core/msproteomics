# nf-core/msproteomics: Database Preparation

## Overview

Protein identification in mass spectrometry proteomics requires a FASTA database of protein sequences to search spectra against.
Every run takes its database from `--database`, as FragPipe headless always searches the workflow's `database.db-path`.
Without `--database`, `--proteome` selects an iProteome from a versioned reference store (see [iProteomes](usage.md#iproteomes)).

## Providing the Database

Provide the FASTA path with `--database`:

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

- **FragPipe mode**: the FASTA is searched as given, so it must already hold decoys whose protein IDs start with the workflow's `database.decoy-tag` (usually `rev_`).
  FragPipe's own decoy check (`FragpipeRun.checkDbConfig`, in FRAGPIPEPLAN's `--write-plan` run, or FragPipe headless) stops the run when the decoys are 0 %, 100 % or outside 40-60 % of the entries.
  To add decoys and contaminants first, set `--add_philosopher_decoy` / `--add_philosopher_contaminants`: the pipeline then runs FragPipe's database update, `philosopher database --custom <fasta> [--nodecoys] [--contam --contamprefix]`, which writes `rev_` decoys.
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
This mapping is handled by `convertOrganismToStandardName()` in `subworkflows/local/utils_nfcore_msproteomics_pipeline/main.nf`.

The following aliases are currently supported:

| Standard Name              | Accepted Aliases                          |
| -------------------------- | ----------------------------------------- |
| _Homo sapiens_             | `human`, `hs`, `Homo_sapiens`             |
| _Mus musculus_             | `mouse`, `mm`, `Mus_musculus`             |
| _Saccharomyces cerevisiae_ | `yeast`, `sc`, `Saccharomyces_cerevisiae` |

For other supported organisms, use the full binomial name (e.g., `Rattus norvegicus`, `Danio rerio`) in your samplesheet.
Underscores and spaces are treated equivalently, and matching is case-insensitive.
