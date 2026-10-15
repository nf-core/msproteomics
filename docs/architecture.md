# nf-core/msproteomics: Workflow Architecture

## Pipeline Overview

nf-core/msproteomics is a Nextflow-based mass spectrometry proteomics pipeline that supports multiple analysis modes.
The entry point (`main.nf`) follows the nf-core template: `PIPELINE_INITIALISATION` (parameter validation against `nextflow_schema.json`, `--help`, `--version`), then `NFCORE_MSPROTEOMICS`, then `PIPELINE_COMPLETION` (summary, email, hooks).
`NFCORE_MSPROTEOMICS` routes execution to the appropriate workflow based on `--mode` and `--tmt_mode`, collects software versions into `pipeline_info/nf_core_msproteomics_software_mqc_versions.yml` and runs MultiQC.
Each mode workflow parses the samplesheet with nf-schema (`assets/schema_input.json`).

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
  ├── GENERATESDRFFROMSAMPLESHEET ─ bookkeeping
  │
  ▼
PARSEFRAGPIPEWORKFLOW ────── read the .workflow file (tool run flags and settings)
  │
  ▼
CHECKCENTROID ────────────── FragPipe's centroid check of the LC-MS files
  │
  ▼
FRAGPIPE_DATABASE ────────── decoy check (CHECKDBCONFIG); optional PHILOSOPHERDATABASE (decoys + contaminants)
  │
  ▼
FRAGPIPE_CONVERT ─────────── optional DIA-Umpire / diaTracer (LC-MS files searched as given)
  │
  ▼
FRAGPIPE_SEARCH ──────────── MSFragger database search (split database: SORTFASTA, MSFRAGGERCALIBRATE, SPLITFASTA, MERGESPLITSEARCH)
  │                           + CrystalC artifact removal (optional)
  │                           + MSBooster rescoring (optional)
  │
  ▼
FRAGPIPE_VALIDATE ────────── PSM validation via Percolator or PeptideProphet
  │                           + PTMProphet for modification localization (optional)
  │
  ▼
FRAGPIPE_INFERENCE ───────── ProteinProphet, Philosopher Filter, iProphet, Abacus
  │
  ▼
FRAGPIPE_PHILOSOPHER_QUANT ─ FreeQuant / LabelQuant, final Philosopher report
  │                           then Export matched fragments (FRAGPIPE_MATCHED_FRAGMENTS)
  │
  ├──► WRITESUBMZML ──────── optional _sub.mzML per LC-MS file; then FRAGPIPESECONDPASS writes the
  │                           second-pass workflow and manifest (workflows/fragpipe.nf)
  ▼
FRAGPIPE_GLYCO ───────────── PairScans, O-Pair, PTM-Shepherd (FRAGPIPE_PTM), MBG
  │
  ▼
FRAGPIPE_QUANT ───────────── IonQuant (LFQ, MS1, isobaric), TMT-Integrator
  │
  ▼
FRAGPIPE_SPECLIB ─────────── spectral library generation (SPECLIBGEN, DIANNFRAGPIPE)
  │
  ▼
FRAGPIPE_EXPORT ──────────── SAINTexpress (optional), FPOP, FPOP coADAPTr, Metaproteomics
```

**FILE_PREPARATION** converts vendor-specific raw files (Thermo `.raw`, Bruker `.d`) to open-format mzML.
Handles decompression and peak picking as needed.

**DIA_PROTEOMICS_ANALYSIS** is an nf-core subworkflow that runs DIA-NN in two stages: first building a spectral library from the FASTA database, then performing the main analysis using that library.
An existing spectral library can be provided via `--diann_speclib` to skip library generation.

**MSstats** performs differential abundance testing when an experimental design with conditions is provided.
Skipped when `--skip_post_msstats` is set.

**IQ** applies the MaxLFQ algorithm to produce protein-level quantification from precursor intensities.

## FragPipe Pipeline Mode

Pipeline mode (`--fragpipe_mode pipeline`) decomposes FragPipe into modular Nextflow subworkflows.
Each subworkflow wraps one or more FragPipe tools and can run independently with full Nextflow caching and resume support.

```
Input CSV
  │
  ├── GENERATESDRFFROMSAMPLESHEET ─ bookkeeping
  │
  ▼
FRAGPIPE_CONVERT ─────────── optional DIA-Umpire / diaTracer (LC-MS files searched as given)
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
  ├──► FRAGPIPE_PHILOSOPHER_QUANT ─ FreeQuant / LabelQuant, final Philosopher report
  │
  ├──► FRAGPIPE_GLYCO ────── PTMShepherd (FRAGPIPE_PTM), O-Pair, MBG glycoproteomics
  │
  ├──► FRAGPIPE_QUANT ────── IonQuant (LFQ) or TMTIntegrator (TMT)
  │
  ├──► FRAGPIPE_SPECLIB ──── Spectral library generation (SPECLIBGEN, DIANNFRAGPIPE)
  │
  └──► FRAGPIPE_EXPORT ───── FPOP, Metaproteomics; SAINTexpress (optional post-run, downstream/)
```

The first step parses the `.workflow` configuration file (`PARSEFRAGPIPEWORKFLOW`) to determine which tools are enabled and extract their parameters.
Only enabled tools actually execute -- the workflow dynamically adapts to whatever the `.workflow` file specifies.

**FRAGPIPE_DATABASE** checks the database decoys the way FragPipe headless does (`CHECKDBCONFIG`); the FASTA is searched as given.
With `--add_philosopher_decoy` / `--add_philosopher_contaminants` it first runs FragPipe's database update (`philosopher database --custom <fasta> [--nodecoys] [--contam --contamprefix]`).

**FRAGPIPE_CONVERT** runs the optional DIA preprocessing (DIA-Umpire, diaTracer); LC-MS files are not converted, because MSFragger reads `.raw`, `.mzML`, `.mzXML`, `.mgf`, `.mzBIN` and Bruker `.d` natively, as in FragPipe.

**FRAGPIPE_SEARCH** runs MSFragger for database searching with native mass calibration.
MSBooster can optionally rescore PSMs using deep learning predictions, and CrystalC removes chimeric spectrum artifacts.

**FRAGPIPE_VALIDATE** performs PSM-level FDR control using either Percolator (semi-supervised learning) or PeptideProphet (mixture modeling).
PTMProphet can optionally localize modification sites.

**FRAGPIPE_INFERENCE** aggregates per-sample results for protein-level inference via ProteinProphet, then applies Philosopher Filter for multi-level FDR control (PSM, peptide, protein).

**Post-inference steps** run in FragPipe's order, each on the tables the previous one left:

- **FRAGPIPE_PHILOSOPHER_QUANT** -- FreeQuant, LabelQuant and the final Philosopher report, then Export matched fragments; Write sub mzML (`tab-run.write_sub_mzml=true`) runs per LC-MS file on the final report
- **FRAGPIPE_GLYCO** -- PAIRSCANS and O-Pair (O-linked), MBG, and `FRAGPIPE_PTM` (PTMShepherd for open search PTM summarization and localization)
- **FRAGPIPE_QUANT** -- IonQuant for label-free quantification (with MaxLFQ and match-between-runs) and isobaric quantification, TMTIntegrator for isobaric quantification
- **FRAGPIPE_SPECLIB** -- spectral library generation for DIA or targeted methods
- **FRAGPIPE_EXPORT** -- export results to FPOP or metaproteomics analysis, and the optional post-run SAINTexpress (FragPipe's Downstream-tab action, written to `downstream/saintexpress/`) (Skyline is not supported: FragPipe 24.0 refuses to run it on Linux, where it finds no Skyline installation)

## FragPipe Headless Mode

Headless mode (`--fragpipe_mode headless`) runs the entire FragPipe pipeline as a single process, exactly reproducing the behavior of `fragpipe --headless` on the command line.

```
Input CSV
  │
  ▼
PARSEFRAGPIPEWORKFLOW ──── read the .workflow file (decoy tag and tool run flags)
  │
  ▼
FRAGPIPE_DATABASE ──────── decoy check (CHECKDBCONFIG); optional PHILOSOPHERDATABASE (decoys + contaminants)
  │
  ▼
FRAGPIPEHEADLESS ───────── all-in-one FragPipe execution on all LC-MS files as given
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
MSPROTEOMICS_FRAGPIPE ───── the label check workflow, as any FragPipe run
  │                         (--fragpipe_mode pipeline: FRAGPIPE; headless: FRAGPIPE_HEADLESS_WF)
  │                         TMT set as VARIABLE modification
  ├── each group's Report psm.tsv ──────── TMTLABELCHECKANALYZE (qc/)
  └── IonQuant combined_modified_peptide ─ TMTLABELCHECKANALYZE_IONQUANT (qc_combined/)
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
  ├── conf/iproteomes.config           ← iProteomes: reference FASTAs and hosted indices (--proteome; skipped with --iproteomes_ignore)
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
