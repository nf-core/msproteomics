# nf-core/msproteomics: Quick Start

Get from raw mass spectrometry data to quantified proteins in minutes.
This guide covers the three main analysis modes: DIA, DDA LFQ, and TMT Label Check.

## Prerequisites

- **Nextflow** >= 25.04.0
- **Docker** or **Singularity**
- **Java** 17+

Install Nextflow if you do not have it:

```bash
curl -s https://get.nextflow.io | bash
```

## Step 1: Install the Pipeline

```bash
nextflow pull nf-core/msproteomics
```

## Step 2: Prepare Your Samplesheet

Create a CSV file with one row per MS run.
The required columns are `sample` and `spectra`; the optional columns are `condition`, `label`, and `fraction`.

| Column      | Required | Description                                          |
| ----------- | -------- | ---------------------------------------------------- |
| `sample`    | Yes      | Unique sample name (no spaces)                       |
| `spectra`   | Yes      | Path to spectra file (`.raw`, `.mzML`, `.d`, `.dia`) |
| `condition` | No       | Experimental group (e.g., control, treated)          |
| `label`     | No       | TMT label channel (e.g., `TMT126`, `TMT127N`)        |
| `fraction`  | No       | Fraction number (positive integer)                   |

Here is a minimal DIA or LFQ samplesheet (`samplesheet.csv`):

```csv
sample,spectra,condition,label,fraction
Sample_A1,/data/raw/A1.raw,control,,
Sample_A2,/data/raw/A2.raw,control,,
Sample_B1,/data/raw/B1.raw,treated,,
Sample_B2,/data/raw/B2.raw,treated,,
```

For a TMT samplesheet example with label annotations, see `assets/samplesheet_tmt.csv`.

## Step 3: Run a DIA Analysis

```bash
nextflow run nf-core/msproteomics \
    --mode diann \
    --input samplesheet.csv \
    --outdir results \
    -profile docker
```

Reference proteomes auto-resolve for common organisms (human, mouse, yeast) based on sample metadata.
For other organisms, supply a FASTA database explicitly with `--database /path/to/database.fasta`.

## Step 4: Run a DDA LFQ Analysis

```bash
nextflow run nf-core/msproteomics \
    --mode fragpipe \
    --input samplesheet.csv \
    --database /path/to/database.fasta \
    --fragpipe_container your-registry/fragpipe:24.0 \
    --fragpipe_workflow LFQ-MBR.workflow \
    --outdir results \
    -profile docker
```

DDA LFQ mode requires a custom FragPipe container with licensed tools.
See the [FragPipe Docker build guide](fragpipe-docker/README.md) for instructions on building the container image.

## Step 5: Run a TMT Label Check

```bash
nextflow run nf-core/msproteomics \
    --mode fragpipe \
    --tmt_mode labelcheck \
    --tmt_type TMT16 \
    --input samplesheet.csv \
    --database /path/to/database.fasta \
    --fragpipe_container your-registry/fragpipe:24.0 \
    --outdir results \
    -profile docker
```

The pipeline auto-selects the correct workflow file based on `--tmt_type`.
Supported types: `TMT6`, `TMT10`, `TMT11`, `TMT16`, `TMT18`, `TMTPRO`.

## Step 6: Check Your Results

All outputs are written under the directory specified by `--outdir`.

**DIA results:**

- `diann/report.stats.tsv` — summary statistics from DIA-NN
- `multiqc/multiqc_report.html` — interactive quality control report

**DDA LFQ results:**

- `ionquant/combined_protein.tsv` — protein-level quantification (MaxLFQ)
- `philosopher_filter/` — filtered PSM and protein reports per sample

**TMT Label Check results:**

- `tmt_labelcheck/tmt_labelcheck_report.html` — labeling efficiency report with per-channel statistics

## Step 7: Resume a Failed Run

If a run fails partway through, fix the issue and resume from the cached results:

```bash
nextflow run nf-core/msproteomics \
    --mode diann \
    --input samplesheet.csv \
    --outdir results \
    -profile docker \
    -resume
```

The `-resume` flag skips tasks that completed successfully and reruns only from the point of failure.

## Next Steps

- [Usage Guide](usage.md) — full parameter reference
- [Output Documentation](output.md) — interpreting results for each workflow
- [Database Preparation](database_guide.md) — building custom FASTA databases
- [Troubleshooting](troubleshooting.md) — common issues and solutions
- [Workflow Architecture](architecture.md) — pipeline design and module structure
