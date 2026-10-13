# Module Documentation

Index of all local module documentation for the nf-core/msproteomics pipeline.

## Data Format Conversion

| Module                                        | Description                                                                      |
| --------------------------------------------- | -------------------------------------------------------------------------------- |
| [thermorawfileparser](thermorawfileparser.md) | Convert Thermo `.raw` files to open mass spectrometry formats (DIA-NN path only) |
| [check_centroid](check_centroid.md)           | FragPipe's centroid check of the first LC-MS file before the search              |
| [tdf2mzml](tdf2mzml.md)                       | Convert Bruker `.d` (TDF) files to mzML format                                   |
| [diatracer](diatracer.md)                     | Convert Bruker `.d` files to mzML for DIA analysis                               |
| [diaumpire](diaumpire.md)                     | Generate pseudo-spectra from DIA mass spectrometry data                          |
| [decompress](decompress.md)                   | Decompress .gz/.tar/.tar.gz/.zip inputs to .mzML, .raw or Bruker `.d`            |
| [mzml_indexing](mzml_indexing.md)             | Index mzML files using OpenMS FileConverter for efficient random access          |
| [openms_peak_picker](openms_peak_picker.md)   | Centroid peak picking using OpenMS PeakPickerHiRes                               |

## Database Preparation

| Module                                          | Description                                                                                                              |
| ----------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| [philosopher_database](philosopher_database.md) | FragPipe database update (decoys + contaminants), only with `--add_philosopher_decoy` / `--add_philosopher_contaminants` |
| [split_fasta](split_fasta.md)                   | Split a FASTA database into chunks for parallel MSFragger search                                                         |

## Database Search

| Module                                        | Description                                                             |
| --------------------------------------------- | ----------------------------------------------------------------------- |
| [msfragger](msfragger.md)                     | Ultrafast peptide identification using MSFragger database search engine |
| [msfragger_calibrate](msfragger_calibrate.md) | Mass calibration and parameter optimization for split-database search   |
| [merge_split_search](merge_split_search.md)   | Merge MSFragger split-database search results per sample                |
| [crystalc](crystalc.md)                       | Remove chimeric artifact PSMs from pepXML search results                |

## PSM Validation

| Module                                                      | Description                                                                  |
| ----------------------------------------------------------- | ---------------------------------------------------------------------------- |
| [msbooster](msbooster.md)                                   | Deep learning-based PSM rescoring using predicted retention time and spectra |
| [percolator](percolator.md)                                 | Semi-supervised PSM validation and FDR control using Percolator              |
| [philosopher_peptideprophet](philosopher_peptideprophet.md) | Statistical validation of peptide assignments using PeptideProphet           |
| [ptmprophet](ptmprophet.md)                                 | PTM site localization using PTMProphet                                       |

## Protein Inference

| Module                                                      | Description                                              |
| ----------------------------------------------------------- | -------------------------------------------------------- |
| [philosopher_proteinprophet](philosopher_proteinprophet.md) | Statistical protein inference using ProteinProphet       |
| [philosopher_filter](philosopher_filter.md)                 | FDR filtering and report generation using Philosopher    |
| [philosopher_report](philosopher_report.md)                 | Generate TSV reports from Philosopher workspace binaries |

## Quantification

| Module                                              | Description                                                                           |
| --------------------------------------------------- | ------------------------------------------------------------------------------------- |
| [ionquant](ionquant.md)                             | Label-free and isobaric quantification with match-between-runs and MaxLFQ             |
| [freequant](freequant.md)                           | Label-free quantification using Philosopher's `freequant` command                     |
| [philosopher_labelquant](philosopher_labelquant.md) | Extract TMT/isobaric reporter ion intensities using Philosopher labelquant            |
| [tmtintegrator](tmtintegrator.md)                   | TMT-Integrator for isobaric labeling quantification                                   |
| [iq](iq.md)                                         | Protein quantification using the iq R package with MaxLFQ algorithm                   |
| [msstats_lfq](msstats_lfq.md)                       | Statistical analysis of label-free quantification data using MSstats                  |
| [diann_fragpipe](diann_fragpipe.md)                 | FragPipe-mode DIA-NN quantification with FragPipe's bundled DIA-NN and its post-steps |

## PTM Analysis

| Module                        | Description                                                              |
| ----------------------------- | ------------------------------------------------------------------------ |
| [ptmshepherd](ptmshepherd.md) | Post-translational modification analysis and mass shift characterization |

## Glycoproteomics

| Module            | Description                                                                  |
| ----------------- | ---------------------------------------------------------------------------- |
| [mbg](mbg.md)     | Mass-based glycoproteomics matching for glycan identification                |
| [opair](opair.md) | O-glycoproteomics analysis for O-glycan identification and site localization |

## QC and Reporting

| Module                            | Description                                                                |
| --------------------------------- | -------------------------------------------------------------------------- |
| [tmtlabelcheck](tmtlabelcheck.md) | Analyze TMT labeling efficiency from PSM data                              |
| [pmultiqc](pmultiqc.md)           | Proteomics-specific QC report generation using the pmultiqc MultiQC plugin |

## Export

| Module                                                  | Description                                                                        |
| ------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| [saintexpress](saintexpress.md)                         | Score protein-protein interactions from AP-MS data using SAINTexpress              |
| [fpop](fpop.md)                                         | Fast Photochemical Oxidation of Proteins (FPOP) analysis                           |
| [export_matched_fragments](export_matched_fragments.md) | FragPipe "Export matched fragments": matched fragment ions per PSM in psm.tsv      |
| [metaproteomics](metaproteomics.md)                     | Metaproteomics database optimization and taxonomy analysis                         |
| [speclibgen](speclibgen.md)                             | Generate spectral libraries with FragPipe's gen_con_spec_lib.py (FragPipe-SpecLib) |

## Workflow Orchestration

| Module                                                | Description                                                     |
| ----------------------------------------------------- | --------------------------------------------------------------- |
| [fragpipe_headless](fragpipe_headless.md)             | Run FragPipe in headless mode as a single all-in-one process    |
| [parse_fragpipe_workflow](parse_fragpipe_workflow.md) | Parse FragPipe workflow files into per-tool configuration files |

## Utility

| Module                            | Description                              |
| --------------------------------- | ---------------------------------------- |
| [generate_sdrf](generate_sdrf.md) | Convert a samplesheet CSV to SDRF format |
