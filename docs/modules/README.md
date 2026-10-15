# Module Documentation

Index of all local module documentation for the nf-core/msproteomics pipeline.
Each page is named after the module directory under `modules/local` (the process name is the upper-cased name).

## Data Format Conversion and Preparation

| Module                                  | Description                                                             |
| --------------------------------------- | ----------------------------------------------------------------------- |
| [checkcentroid](checkcentroid.md)       | FragPipe's centroid check of the first LC-MS file before the search     |
| [tdf2mzml](tdf2mzml.md)                 | Convert Bruker `.d` (TDF) files to mzML format                          |
| [diatracer](diatracer.md)               | diaTracer: Bruker diaPASEF `.d` files to mzML for DIA analysis          |
| [diaumpire](diaumpire.md)               | Generate pseudo-spectra from DIA mass spectrometry data                 |
| [decompress](decompress.md)             | Decompress .gz/.tar/.tar.gz/.zip inputs to .mzML, .raw or Bruker `.d`   |
| [mzmlindexing](mzmlindexing.md)         | Index mzML files using OpenMS FileConverter for efficient random access |
| [openmspeakpicker](openmspeakpicker.md) | Centroid peak picking using OpenMS PeakPickerHiRes                      |

## Database Preparation

| Module                                        | Description                                                                                                              |
| --------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| [checkdbconfig](checkdbconfig.md)             | Check the database decoys the way FragPipe headless does                                                                 |
| [philosopherdatabase](philosopherdatabase.md) | FragPipe database update (decoys + contaminants), only with `--add_philosopher_decoy` / `--add_philosopher_contaminants` |
| [splitfasta](splitfasta.md)                   | Split a FASTA database into chunks for parallel MSFragger search                                                         |
| [sortfasta](sortfasta.md)                     | Sort the proteins of a FASTA file for the calibration search of a split-database run                                     |

## Database Search

| Module                                                    | Description                                                                                      |
| --------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| [msfragger](msfragger.md)                                 | Ultrafast peptide identification using MSFragger database search engine                          |
| [msfraggercalibrate](msfraggercalibrate.md)               | Mass calibration and parameter optimization for split-database search                            |
| [msfraggerapplycalibration](msfraggerapplycalibration.md) | Apply MSFragger's calibration output to the search params and spectra of a split-database search |
| [mergesplitsearch](mergesplitsearch.md)                   | Merge MSFragger split-database search results per sample                                         |
| [crystalc](crystalc.md)                                   | Remove chimeric artifact PSMs from pepXML search results                                         |

## PSM Validation

| Module                                                    | Description                                                                  |
| --------------------------------------------------------- | ---------------------------------------------------------------------------- |
| [msbooster](msbooster.md)                                 | Deep learning-based PSM rescoring using predicted retention time and spectra |
| [percolator](percolator.md)                               | Semi-supervised PSM validation and FDR control using Percolator              |
| [philosopherpeptideprophet](philosopherpeptideprophet.md) | Statistical validation of peptide assignments using PeptideProphet           |
| [ptmprophet](ptmprophet.md)                               | PTM site localization using PTMProphet                                       |

## Protein Inference

| Module                                                    | Description                                              |
| --------------------------------------------------------- | -------------------------------------------------------- |
| [philosopheriprophet](philosopheriprophet.md)             | Philosopher iProphet over all validated pepXMLs          |
| [philosopherproteinprophet](philosopherproteinprophet.md) | Statistical protein inference using ProteinProphet       |
| [philosopherfilter](philosopherfilter.md)                 | FDR filtering and report generation using Philosopher    |
| [philosopherreport](philosopherreport.md)                 | Generate TSV reports from Philosopher workspace binaries |
| [philosopherabacus](philosopherabacus.md)                 | Philosopher abacus multi-experiment report               |

## Quantification

| Module                                            | Description                                                                           |
| ------------------------------------------------- | ------------------------------------------------------------------------------------- |
| [ionquant](ionquant.md)                           | Label-free and isobaric quantification with match-between-runs and MaxLFQ             |
| [freequant](freequant.md)                         | Label-free quantification using Philosopher's `freequant` command                     |
| [philosopherlabelquant](philosopherlabelquant.md) | Extract TMT/isobaric reporter ion intensities using Philosopher labelquant            |
| [tmtintegrator](tmtintegrator.md)                 | TMT-Integrator for isobaric labeling quantification                                   |
| [iq](iq.md)                                       | Protein quantification using the iq R package with MaxLFQ algorithm                   |
| [msstatslfq](msstatslfq.md)                       | Statistical analysis of label-free quantification data using MSstats                  |
| [diannfragpipe](diannfragpipe.md)                 | FragPipe-mode DIA-NN quantification with FragPipe's bundled DIA-NN and its post-steps |

## PTM Analysis

| Module                        | Description                                                              |
| ----------------------------- | ------------------------------------------------------------------------ |
| [ptmshepherd](ptmshepherd.md) | Post-translational modification analysis and mass shift characterization |

## Glycoproteomics

| Module                    | Description                                                                  |
| ------------------------- | ---------------------------------------------------------------------------- |
| [pairscans](pairscans.md) | Pair parent and child scans of one LC-MS file for O-Pair                     |
| [mbg](mbg.md)             | Mass-based glycoproteomics matching for glycan identification                |
| [opair](opair.md)         | O-glycoproteomics analysis for O-glycan identification and site localization |

## QC and Reporting

| Module                                          | Description                                                                |
| ----------------------------------------------- | -------------------------------------------------------------------------- |
| [tmtlabelcheckanalyze](tmtlabelcheckanalyze.md) | Analyze TMT labeling efficiency from PSM data                              |
| [pmultiqc](pmultiqc.md)                         | Proteomics-specific QC report generation using the pmultiqc MultiQC plugin |

## Export

| Module                                              | Description                                                                        |
| --------------------------------------------------- | ---------------------------------------------------------------------------------- |
| [saintexpress](saintexpress.md)                     | Score protein-protein interactions from AP-MS data using SAINTexpress              |
| [reprinttosaintinput](reprinttosaintinput.md)       | SAINTexpress inputs from a REPRINT table                                           |
| [fpop](fpop.md)                                     | Fast Photochemical Oxidation of Proteins (FPOP) analysis                           |
| [fpopcoadaptr](fpopcoadaptr.md)                     | FPOP coADAPTr converter over the whole FragPipe work dir                           |
| [exportmatchedfragments](exportmatchedfragments.md) | FragPipe "Export matched fragments": matched fragment ions per PSM in psm.tsv      |
| [metaproteomics](metaproteomics.md)                 | Metaproteomics database optimization and taxonomy analysis                         |
| [speclibgen](speclibgen.md)                         | Generate spectral libraries with FragPipe's gen_con_spec_lib.py (FragPipe-SpecLib) |

## Workflow Orchestration

| Module                                            | Description                                                             |
| ------------------------------------------------- | ----------------------------------------------------------------------- |
| [fragpipeheadless](fragpipeheadless.md)           | Run FragPipe in headless mode as a single all-in-one process            |
| [parsefragpipeworkflow](parsefragpipeworkflow.md) | Parse FragPipe workflow files into per-tool configuration files         |
| [fragpipesaveworkflow](fragpipesaveworkflow.md)   | Write FragPipe's saved UI state (`fragpipe.workflow`) for a run         |
| [fragpipesecondpass](fragpipesecondpass.md)       | FragPipe's second-pass workflow and manifest after a Write sub mzML run |
| [writesubmzml](writesubmzml.md)                   | FragPipe "Write sub mzML": the scans no confident PSM explains          |

## Utility

| Module                                                        | Description                                                     |
| ------------------------------------------------------------- | --------------------------------------------------------------- |
| [generatesdrffromsamplesheet](generatesdrffromsamplesheet.md) | Convert a samplesheet CSV to SDRF format                        |
| [fragpipesdrf](fragpipesdrf.md)                               | Write the FragPipe SDRF template (`sdrf.tsv`) for a run         |
| [writetextfile](writetextfile.md)                             | Write a text file whose body was built outside the tool modules |
