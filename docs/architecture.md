# nf-core/msproteomics: Workflow Architecture

## Pipeline Overview

nf-core/msproteomics is a Nextflow-based mass spectrometry proteomics pipeline that supports multiple analysis modes.
The entry point (`main.nf`) routes execution to the appropriate workflow based on two parameters: `--mode` and `--tmt_mode`.

```
--mode diann                              → DIA Workflow (DIA-NN)
--mode fragpipe                           → FragPipe Workflow (DDA LFQ or TMT Quant)
--mode fragpipe --tmt_mode labelcheck     → TMT Label Check Workflow
--mode fragpipe --tmt_mode quant          → TMT Quantification Workflow
```

Within the FragPipe workflows, the `--fragpipe_mode` parameter controls how tools are executed:

- `pipeline` (default) -- modular subworkflows orchestrated by Nextflow, with per-tool caching and resume
- `headless` -- FragPipe runs as a single all-in-one process, reproducing exact FragPipe CLI behavior

## DIA Workflow

The DIA workflow (`workflows/diann.nf`) runs DIA-NN-based analysis with optional statistical post-processing.

```
Input CSV
  │
  ├── GENERATE_SDRF ─────────── bookkeeping (SDRF metadata generation)
  │
  ▼
FILE_PREPARATION ──────────── RAW → mzML conversion, decompression, peak picking
  │
  ▼
DIA_PROTEOMICS_ANALYSIS ───── DIA-NN library generation + quantification
  │                            (nf-core subworkflow)
  │
  ├──► MSstats (optional) ──── differential abundance analysis
  │
  ├──► DIANNR ─────────────── R-based post-processing and QC
  │
  └──► IQ / MaxLFQ ────────── protein-level quantification via MaxLFQ algorithm
```

**FILE_PREPARATION** converts vendor-specific raw files (Thermo `.raw`, Bruker `.d`) to open-format mzML.
Handles decompression and peak picking as needed.

**DIA_PROTEOMICS_ANALYSIS** is an nf-core subworkflow that runs DIA-NN in two stages: first building a spectral library from the FASTA database, then performing the main analysis using that library.
An existing spectral library can be provided via `--diann_speclib` to skip library generation.

**MSstats** performs differential abundance testing when an experimental design with conditions is provided.
Skipped when `--skip_post_msstats` is set.

**DIANNR** runs R-based post-processing on the DIA-NN report, applying quality filters and generating summary tables.

**IQ** applies the MaxLFQ algorithm to produce protein-level quantification from precursor intensities.

## FragPipe Pipeline Mode

Pipeline mode (`--fragpipe_mode pipeline`) decomposes FragPipe into modular Nextflow subworkflows.
Each subworkflow wraps one or more FragPipe tools and can run independently with full Nextflow caching and resume support.

```
Input CSV
  │
  ├── GENERATE_SDRF ──────── bookkeeping
  │
  ▼
FRAGPIPE_CONVERT ─────────── RAW → mzML (ThermoRawFileParser)
  │                           Optional: DIA-Umpire / DIA-Tracer preprocessing
  │
  │     FRAGPIPE_DATABASE ─── Philosopher: add decoys + contaminants to FASTA
  │         │
  ▼         ▼
FRAGPIPE_SEARCH ──────────── MSFragger database search
  │                           + MSBooster AI rescoring (optional)
  │                           + CrystalC artifact removal (optional)
  │
  ▼
FRAGPIPE_VALIDATE ────────── PSM validation via Percolator or PeptideProphet
  │                           + PTMProphet for modification localization (optional)
  │
  ▼
FRAGPIPE_INFERENCE ───────── ProteinProphet protein-level inference
  │                           + Philosopher Filter FDR filtering
  │
  ├──► FRAGPIPE_QUANT ────── IonQuant (LFQ) or TMTIntegrator (TMT) or FreeQuant
  │
  ├──► FRAGPIPE_PTM ──────── PTMShepherd post-translational modification analysis
  │
  ├──► FRAGPIPE_SPECLIB ──── Spectral library generation (EasyPQP / SpecLibGen)
  │
  ├──► FRAGPIPE_GLYCO ────── Glycoproteomics: MSFragger-Glyco / O-Pair
  │
  └──► FRAGPIPE_EXPORT ───── Skyline, SAINTexpress, FPOP, Metaproteomics
```

The first step parses the `.workflow` configuration file (`PARSE_FRAGPIPE_WORKFLOW`) to determine which tools are enabled and extract their parameters.
Only enabled tools actually execute -- the workflow dynamically adapts to whatever the `.workflow` file specifies.

**FRAGPIPE_DATABASE** runs Philosopher to append decoy sequences and common contaminants to the protein FASTA database.
If a prebuilt database is provided, this step is skipped.

**FRAGPIPE_CONVERT** converts Thermo `.raw` files to mzML using ThermoRawFileParser.
For DIA data, optional DIA-Umpire or DIA-Tracer preprocessing runs after conversion.

**FRAGPIPE_SEARCH** runs MSFragger for database searching with native mass calibration.
MSBooster can optionally rescore PSMs using deep learning predictions, and CrystalC removes chimeric spectrum artifacts.

**FRAGPIPE_VALIDATE** performs PSM-level FDR control using either Percolator (semi-supervised learning) or PeptideProphet (mixture modeling).
PTMProphet can optionally localize modification sites.

**FRAGPIPE_INFERENCE** aggregates per-sample results for protein-level inference via ProteinProphet, then applies Philosopher Filter for multi-level FDR control (PSM, peptide, protein).

**Post-inference steps** run in parallel after inference completes:

- **FRAGPIPE_QUANT** -- IonQuant for label-free quantification (with MaxLFQ and match-between-runs), TMTIntegrator for isobaric quantification, or FreeQuant for spectral counting
- **FRAGPIPE_PTM** -- PTMShepherd for open search PTM summarization and localization
- **FRAGPIPE_SPECLIB** -- spectral library generation for DIA or targeted methods
- **FRAGPIPE_GLYCO** -- MSFragger-Glyco (N-linked) and O-Pair (O-linked) glycoproteomics
- **FRAGPIPE_EXPORT** -- export results to Skyline, SAINTexpress, or other downstream tools

## FragPipe Headless Mode

Headless mode (`--fragpipe_mode headless`) runs the entire FragPipe pipeline as a single process, exactly reproducing the behavior of `fragpipe --headless` on the command line.

```
Input CSV
  │
  ▼
THERMORAWFILEPARSER ────── RAW → mzML conversion (per-sample, parallel)
  │
  │   PHILOSOPHER_DATABASE ── decoy generation (required for headless)
  │       │
  ▼       ▼
FRAGPIPE_HEADLESS ─────── all-in-one FragPipe execution
                            (single process, all tools run sequentially)
```

Headless mode is useful when:

- The `.workflow` file uses tool combinations not yet supported by pipeline mode
- You need to reproduce exact FragPipe CLI behavior for validation
- Debugging differences between pipeline mode and standalone FragPipe

The trade-off is that headless mode runs as a single process, so there is no per-tool Nextflow caching or resume.
If the process fails partway through, it must restart from the beginning.

A manifest file is auto-generated from the samplesheet to organize samples into experiments.
For TMT experiments, per-plex annotation files are generated automatically from the `label` column in the samplesheet.

## TMT Label Check Workflow

The TMT Label Check workflow (`workflows/tmt_labelcheck.nf`) is a QC workflow that measures TMT labeling efficiency before committing to a full quantification experiment.

```
Input CSV
  │
  ├── GENERATE_SDRF ──────── bookkeeping
  │
  ▼
FRAGPIPE_DATABASE ────────── prepare FASTA with decoys
  │
  ▼
FRAGPIPE_CONVERT ─────────── RAW → mzML conversion
  │
  ▼
MSFRAGGER ────────────────── database search (single-pass, no calibration)
  │                           TMT set as VARIABLE modification
  │
  ▼
FRAGPIPE_VALIDATE ────────── Percolator PSM-level FDR
  │
  ▼
FRAGPIPE_INFERENCE ───────── ProteinProphet + Philosopher Filter
  │
  ▼
FRAGPIPE_QUANT ───────────── IonQuant (aggregate PSMs)
  │
  ▼
TMT_LABELCHECK_ANALYZE ───── labeling efficiency report (HTML + TSV)
```

The key difference from a standard search is that TMT is configured as a **variable** modification (not fixed).
This allows detection of both labeled and unlabeled peptides.
Labeling efficiency is calculated as:

```
Efficiency = sum(Labeled_Sites) / sum(Total_Sites) x 100
```

Where `Total_Sites` counts all lysines plus the peptide N-terminus, and `Labeled_Sites` counts sites carrying a TMT modification.

The workflow auto-selects the appropriate `.workflow` file from `assets/TMT-labelcheck-*.workflow` based on `--tmt_type` (TMT6, TMT10, TMT11, TMT16, TMT18, or TMTPRO).

## TMT Quantification Workflow

TMT quantification uses the same FragPipe pipeline mode described above, with TMTIntegrator enabled in the `.workflow` file.

When samples have a `label` column in the samplesheet, the workflow automatically:

1. Builds a TMT annotation file mapping each channel to its sample name, condition, and replicate
2. Passes the annotation to `FRAGPIPE_QUANT`, which runs TMTIntegrator for isobaric quantification
3. Produces normalized abundance tables at the protein, peptide, and gene levels

When no `label` column is present, the annotation is empty and FragPipe runs in LFQ mode instead.

## Configuration Loading

Configurations are loaded in a specific order, with later files overriding earlier ones.

```
nextflow.config                        ← root params, profiles, plugins (always loaded)
  │
  ├── conf/base.config                 ← resource labels (process_low, process_high, etc.)
  │
  ├── conf/modules.config              ← per-module ext.args, publishDir, resource overrides
  │
  ├── conf/reference_proteomes.config  ← pre-built organism databases (human, mouse, yeast)
  │
  ├── conf/instruments/*.config        ← instrument-specific overlays (loaded via -c flag)
  │     e.g., thermo_ascend, bruker_flex
  │
  └── conf/variants/*.config           ← analysis variant overlays (loaded via -c flag)
        e.g., diaphos (DIA phosphoproteomics)
```

**Resource labels** (`conf/base.config`) define standard resource tiers (`process_low`, `process_medium`, `process_high`, etc.) that modules reference.
Override these in your profile or custom config to match your compute environment.

**Module options** (`conf/modules.config`) set per-module CLI arguments via `ext.args`, output directories via `publishDir`, and resource overrides.

**Instrument overlays** adjust mass accuracy, m/z ranges, and scan window parameters for specific mass spectrometers.
Load them with `-c conf/instruments/thermo_ascend.config`.

**Variant overlays** configure specialized analysis modes such as phosphoproteomics.
Load them with `-c conf/variants/diaphos.config`.

## Adding a New Instrument

To add support for a new mass spectrometer:

1. Copy an existing instrument config (e.g., `conf/instruments/thermo_ascend.config`)
2. Adjust parameters for your instrument:
   - Mass accuracy (`--mass-acc` for DIA-NN, `precursor_mass_lower/upper` for MSFragger)
   - m/z scan range
   - Scan window settings
   - Ion mobility parameters (for trapped ion mobility instruments)
3. Save as `conf/instruments/<vendor>_<model>.config`
4. Load at runtime: `nextflow run main.nf -c conf/instruments/<vendor>_<model>.config ...`

See [docs/usage.md](usage.md) for the full parameter reference table.
