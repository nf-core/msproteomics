# nf-core/msproteomics: Usage

## Introduction

nf-core/msproteomics is a mass spectrometry proteomics preprocessing pipeline supporting DIA, DDA LFQ, TMT label check, and generic FragPipe workflows.
This document describes how to prepare input data, configure the pipeline, and run each workflow type.

## Input Format

### CSV Samplesheet

The pipeline accepts a CSV samplesheet describing your samples and their raw data file locations.

```bash
--input samplesheet.csv
```

The samplesheet must contain the following columns:

| Column | Description | Required |
| --- | --- | --- |
| `sample` | Unique sample identifier | Yes |
| `spectra` | Path or URI to the raw/mzML data file | Yes |
| `condition` | Experimental condition or group (used for downstream statistics; defaults to sample name) | No |
| `label` | Isobaric label channel (e.g., TMT126, TMT127N); leave empty for LFQ/DIA | No |
| `fraction` | Fraction number for fractionated experiments; leave empty for single-shot | No |

Example samplesheet (DIA or DDA LFQ):

```csv
sample,spectra,condition,label,fraction
Sample_A1,/data/raw/A1.raw,control,,
Sample_A2,/data/raw/A2.raw,control,,
Sample_B1,/data/raw/B1.raw,treated,,
Sample_B2,/data/raw/B2.raw,treated,,
```

Example samplesheet (TMT):

```csv
sample,spectra,condition,label,fraction
Sample_1,/data/raw/pool1.raw,control,TMT126,
Sample_1,/data/raw/pool1.raw,control,TMT127N,
Sample_1,/data/raw/pool1.raw,treated,TMT128C,
```

The pipeline internally generates an SDRF file from the samplesheet via the `GENERATE_SDRF_FROM_SAMPLESHEET` module.
This SDRF is used for bookkeeping and compatibility with downstream tools that expect SDRF-formatted metadata.

## Analysis Modes

The `--mode` parameter selects the analysis engine:

| Mode | Description | Engine |
| --- | --- | --- |
| `diann` | DIA quantitative proteomics | [DIA-NN](https://github.com/vdemichev/DiaNN) |
| `fragpipe` | DDA LFQ, TMT, or generic FragPipe workflows | [FragPipe](https://fragpipe.nesvilab.org/) |

### DIA Mode (`--mode diann`)

DIA workflows use DIA-NN for library-free or library-based DIA analysis.
DIA-specific variants (phospho) are applied via `-c conf/variants/*.config` files (see [DIA Method Variants](#dia-method-variants)).

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -profile docker
```

Reference proteomes can be automatically fetched from UniProt for supported organisms when configured via `conf/reference_proteomes.config`.
To use a custom database, provide `--database /path/to/database.fasta`.

### FragPipe Mode (`--mode fragpipe`)

FragPipe mode supports DDA LFQ, TMT label check, TMT quantification, and generic FragPipe workflows.
The specific workflow is controlled by the `--fragpipe_workflow` and `--tmt_mode` parameters.

**DDA LFQ** -- uses MSFragger, MSBooster, Percolator, ProteinProphet, Philosopher, and IonQuant:

```bash
nextflow run nf-core/msproteomics \
  --mode fragpipe \
  --fragpipe_workflow LFQ-MBR.workflow \
  --input samplesheet.csv \
  --database /path/to/uniprot_human.fasta \
  --outdir results \
  -profile docker
```

**TMT Label Check** -- assess TMT labeling efficiency:

```bash
nextflow run nf-core/msproteomics \
  --mode fragpipe \
  --tmt_mode labelcheck \
  --input samplesheet.csv \
  --database /path/to/uniprot_human.fasta \
  --outdir results \
  -profile docker
```

**TMT Quantification** -- full TMT isobaric quantification:

```bash
nextflow run nf-core/msproteomics \
  --mode fragpipe \
  --tmt_mode quant \
  --input samplesheet.csv \
  --database /path/to/uniprot_human.fasta \
  --outdir results \
  -profile docker
```

**Generic FragPipe** -- run any FragPipe workflow by providing a `.workflow` configuration file:

```bash
nextflow run nf-core/msproteomics \
  --mode fragpipe \
  --fragpipe_workflow /path/to/custom.workflow \
  --input samplesheet.csv \
  --database /path/to/uniprot_human.fasta \
  --outdir results \
  -profile docker
```

The `--fragpipe_mode` parameter controls execution mode:

- `pipeline` (default): Modular Nextflow subworkflow execution
- `allinone`: Single-process FragPipe headless execution

## Key Parameters

| Parameter | Description | Default |
| --- | --- | --- |
| `--input` | Path to CSV samplesheet | Required |
| `--mode` | Analysis mode: `diann` or `fragpipe` | Required |
| `--outdir` | Output directory | Required |
| `--database` | FASTA protein database | Auto-selected for DIA; required for FragPipe workflows |
| `--fragpipe_workflow` | FragPipe `.workflow` file or workflow name (FragPipe mode) | None |
| `--tmt_mode` | TMT analysis mode: `labelcheck` or `quant` (FragPipe mode) | None |
| `--fragpipe_mode` | FragPipe execution mode: `pipeline` or `allinone` | `pipeline` |

## Container Requirements

### DIA Workflows

DIA workflows use publicly available containers for DIA-NN and associated tools.
No special licensing is required.

### FragPipe-Based Workflows (DDA LFQ, TMT, Generic)

The pipeline defaults to the public base image `docker.io/fcyucn/fragpipe:24.0`, which includes open-source FragPipe tools (Philosopher, Percolator, MSBooster, PTMShepherd, etc.).

However, commercial/licensed tools (MSFragger, IonQuant, DiaTracer) are not included in the public image.
You must build your own container image with these tools under your own license terms.
See [Docker build instructions](fragpipe-docker/README.md) for academic and commercial container builds.

**FragPipe is free for academic use.**
Commercial users must obtain a license from Fragmatics.

## Profiles

Use `-profile` to select a software packaging method:

| Profile | Description |
| --- | --- |
| `docker` | Run with Docker containers (recommended) |
| `singularity` | Run with Singularity containers |
| `apptainer` | Run with Apptainer containers |
| `podman` | Run with Podman containers |
| `conda` | Run with Conda environments (last resort) |
| `test` | Minimal stub test for CI validation (DIA mode) |
| `test_dia` | DIA workflow test with small DIA dataset |
| `test_dda_lfq` | DDA LFQ workflow test with small DDA dataset |
| `test_tmt` | TMT label check workflow test |
| `test_tmtq` | TMT quantification workflow test |

Multiple profiles can be combined: `-profile test_dia,docker`

## Running the Pipeline

### Basic Execution

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -profile docker
```

### Using a Parameters File

```bash
nextflow run nf-core/msproteomics -profile docker -params-file params.yaml
```

```yaml
input: 'samplesheet.csv'
mode: 'fragpipe'
fragpipe_workflow: 'LFQ-MBR.workflow'
database: '/path/to/uniprot_human.fasta'
outdir: 'results'
```

### Resuming a Failed Run

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -profile docker \
  -resume
```

### Updating the Pipeline

```bash
nextflow pull nf-core/msproteomics
```

### Reproducibility

Specify a pipeline version with `-r` to ensure reproducible results:

```bash
nextflow run nf-core/msproteomics -r 1.0.0 \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -profile docker
```

## Configuration Hierarchy

The pipeline loads configuration files automatically based on the `--mode` parameter:

**DIA workflows** (`--mode diann`):

1. `nextflow.config` -- root config, always loaded
2. `conf/base_configs/diann.config` -- DIA-NN computed params and helper functions
3. `conf/variants/*.config` -- variant-specific overrides (diaphos) applied via `-c`
4. `conf/instruments/*.config` -- instrument-specific overrides applied via `-c`

**FragPipe workflows** (`--mode fragpipe`):

1. `nextflow.config` -- root config, always loaded (all params defined here)
2. `--tmt_type` selects the TMT plex (TMT6, TMT10, TMT11, TMT16, TMT18, TMTPRO); the correct labelcheck workflow file is auto-selected
3. `conf/instruments/*.config` -- instrument-specific overrides applied via `-c`

## Custom Configuration

### Resource Requests

To change compute resources for specific processes, see the [nf-core resource tuning documentation](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources).

### Custom Tool Arguments

To pass additional arguments to specific tools, use the `ext.args` directive in a custom config file.
See the [nf-core tool arguments documentation](https://nf-co.re/docs/usage/configuration#customising-tool-arguments).

> [!WARNING]
> Do not use `-c <file>` to specify pipeline parameters.
> Custom config files specified with `-c` must only be used for resource tuning, output directories, or module arguments (`ext.args`).

### Institutional Configs

The pipeline dynamically loads configurations from [nf-core/configs](https://github.com/nf-core/configs).
Check if your institution has a pre-configured profile available.

## Running on Seqera Platform (Tower)

### Nextflow Version

Set the Nextflow version in your launch pre-run script:

```bash
export NXF_VER=25.10.2
```

### Compute Environments

FragPipe processes require compute environments with sufficient memory.
Recommended instance families: `r6i` (memory-optimized).

## Instrument-Specific Settings

DIA-NN in-silico library generation uses different m/z ranges and mass accuracy settings depending on the mass spectrometer.
The pipeline ships per-instrument config files in `conf/instruments/` that override the defaults when included with `-c`.

| Abbreviation | Instrument | Config File | diann_min_pr_mz | diann_max_pr_mz | diann_min_fr_mz | diann_max_fr_mz | diann_library_mass_acc | diann_library_ms1_acc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ASC | Thermo Ascend | `conf/instruments/thermo_ascend.config` (default) | 350 | 1050 | 200 | 1800 | 18 | 5 |
| FLX | Bruker timsTOF Flex | `conf/instruments/bruker_flex.config` | 100 | 1700 | 100 | 1700 | 15 | 15 |

Usage example:

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -c conf/instruments/bruker_flex.config \
  -profile docker
```

Thermo Ascend users can omit the `-c` flag since the defaults are already optimized for Ascend instruments.

## DIA Method Variants

Specialized DIA analysis modes are available as variant config files in `conf/variants/`.
These configs adjust DIA-NN parameters (e.g., variable modifications, phosphosite monitoring) for specific experimental designs.

| Variant | Config File | Description |
| --- | --- | --- |
| DIA Phosphoproteomics | `conf/variants/diaphos.config` | Phosphosite localization with UniMod:21 monitoring |

Usage example:

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -c conf/variants/diaphos.config \
  -profile docker
```

Variant configs can be combined with instrument configs by specifying multiple `-c` flags:

```bash
nextflow run nf-core/msproteomics \
  --mode diann \
  --input samplesheet.csv \
  --database /path/to/database.fasta \
  --outdir results \
  -c conf/instruments/bruker_flex.config \
  -c conf/variants/diaphos.config \
  -profile docker
```

## Nextflow Memory Requirements

To limit Nextflow's JVM memory usage, add the following to your environment:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
