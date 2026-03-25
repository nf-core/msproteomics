/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UNIFIED FRAGPIPE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A dynamic nf-core compliant subworkflow that runs FragPipe tools based on
    workflow configuration.

    This subworkflow composes focused subworkflows:
    1. PARSE_FRAGPIPE_WORKFLOW - Parse workflow file for tool configs
    2. FRAGPIPE_DATABASE - Database preparation (or use prebuilt)
    3. FRAGPIPE_CONVERT - File format conversion
    4. FRAGPIPE_SEARCH - MSFragger search with native calibration
    5. FRAGPIPE_VALIDATE - Percolator / PeptideProphet / PTMProphet
    6. FRAGPIPE_INFERENCE - ProteinProphet + Philosopher Filter
    7. FRAGPIPE_PTM - PTMShepherd (optional)
    8. FRAGPIPE_QUANT - IonQuant / TMTIntegrator / FreeQuant
    9. FRAGPIPE_SPECLIB - Spectral library generation (optional)
    10. FRAGPIPE_GLYCO - Glycoproteomics (optional)
    11. FRAGPIPE_EXPORT - Skyline / SAINT / FPOP / Metaproteomics (optional)

    Tool Execution Order (DAG):
    DATABASE + CONVERT (parallel)
        -> SEARCH (MSFragger + optional CrystalC/MSBooster)
        -> VALIDATE (Percolator or PeptideProphet, optional PTMProphet)
        -> INFERENCE (ProteinProphet + Philosopher Filter)
        -> PTM / QUANT / SPECLIB / GLYCO / EXPORT (parallel post-inference)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import config parser
include { PARSE_FRAGPIPE_WORKFLOW  } from '../../../modules/local/parse_fragpipe_workflow/main'

// Import utility functions
include { shouldRunTool } from '../fragpipe_utils'

// Import subworkflows
include { FRAGPIPE_DATABASE                                    } from '../fragpipe_database/main'
include { FRAGPIPE_CONVERT                                     } from '../fragpipe_convert/main'
include { FRAGPIPE_SEARCH                                      } from '../fragpipe_search/main'
include { FRAGPIPE_VALIDATE  } from '../fragpipe_validate/main'
include { FRAGPIPE_INFERENCE } from '../fragpipe_inference/main'
include { FRAGPIPE_PTM       } from '../fragpipe_ptm/main'
include { FRAGPIPE_QUANT     } from '../fragpipe_quant/main'
include { FRAGPIPE_SPECLIB   } from '../fragpipe_speclib/main'
include { FRAGPIPE_GLYCO     } from '../fragpipe_glyco/main'
include { FRAGPIPE_EXPORT    } from '../fragpipe_export/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FRAGPIPE_WF {
    take:
    ch_file            // channel: [ val(sample_meta), path(mzml/raw/d) ] - sample files
    ch_database        // channel: [ val(meta), path(fasta) ] - protein database (empty if prebuilt)
    ch_workflow_file   // channel: [ val(meta), path(workflow_file) ] - FragPipe workflow file
    ch_decoy_db        // channel: [ val(meta), path(fasta) ] - prebuilt philosopher database (empty if not provided)
    ch_annotation      // channel: path(annotation_file) - value channel, [] when not used
    ch_manifest        // channel: path(manifest.fp-manifest) - value channel, [] when not used

    main:

    //
    // STEP 1: Parse workflow config to get tool run flags and args
    //
    PARSE_FRAGPIPE_WORKFLOW(ch_workflow_file)

    // Read tool configs JSON as a string value channel.
    // Strings avoid the Groovy Map-as-named-parameter bug in channel operators.
    // Use splitText to avoid blocking I/O inside map{} on cloud executors.
    ch_tool_configs = PARSE_FRAGPIPE_WORKFLOW.out.tool_configs_json
        .map { _meta, json_file -> json_file }
        .splitText()
        .collect()
        .map { lines -> lines.join('') }
        .first()

    // Extract config files from PARSE_FRAGPIPE_WORKFLOW as value channels.
    // These are params files for tools that use them (MSFragger, MSBooster, etc.).
    // The parse module touches fragger.params, msbooster.config, peptide-prophet.config
    // so those always exist. tmtintegrator.yml is truly optional.
    ch_fragger_params = PARSE_FRAGPIPE_WORKFLOW.out.msfragger_config
        .map { _meta, params -> params }
        .first()

    ch_msbooster_params = PARSE_FRAGPIPE_WORKFLOW.out.msbooster_config
        .map { _meta, params -> params }
        .first()

    ch_peptideprophet_config = PARSE_FRAGPIPE_WORKFLOW.out.peptideprophet_config
        .map { _meta, config -> config }
        .first()

    ch_tmtintegrator_config = PARSE_FRAGPIPE_WORKFLOW.out.tmtintegrator_config
        .map { _meta, config -> config }
        .first()

    // Determine data type from workflow config (DIA if diaumpire/diatracer enabled, else DDA)
    // ch_tool_configs is a value channel (via .first()), so .map always emits exactly one item.
    ch_data_type = ch_tool_configs
        .map { configs ->
            (shouldRunTool(configs, 'diaumpire') || shouldRunTool(configs, 'diatracer')) ? 'DIA' : 'DDA'
        }

    //
    // STEP 2: Prepare database with decoys and contaminants
    // If prebuilt database (ch_decoy_db) is provided, use it directly.
    // Otherwise, run PHILOSOPHER_DATABASE to add decoys and contaminants.
    // concat().first(): prebuilt items come first; if empty, falls through to generated.
    //
    FRAGPIPE_DATABASE(ch_database)
    ch_fasta = ch_decoy_db
        .map { _meta, fasta -> fasta }
        .concat(FRAGPIPE_DATABASE.out.fasta_path)
        .first()

    //
    // STEP 3: Convert RAW to mzML + optional DIA preprocessing (per-sample, parallel)
    //
    FRAGPIPE_CONVERT(ch_file, ch_tool_configs)

    ch_mzml = FRAGPIPE_CONVERT.out.mzml

    //
    // STEP 4: Database search (MSFragger with native calibration)
    //
    FRAGPIPE_SEARCH(
        ch_mzml,
        ch_fasta,
        ch_tool_configs,
        ch_fragger_params,
        ch_msbooster_params,
        params.msfragger_num_slices
    )
    ch_search_pepxml       = FRAGPIPE_SEARCH.out.pepxml
    ch_search_pin          = FRAGPIPE_SEARCH.out.pin
    ch_search_pin_rescored = FRAGPIPE_SEARCH.out.pin_rescored

    // Original spectra for downstream tools (IonQuant needs MS1 scans).
    ch_mzml_downstream = ch_mzml

    //
    // STEP 5: Percolator / PeptideProphet / PTMProphet validation (per-sample, parallel)
    //
    FRAGPIPE_VALIDATE(
        ch_search_pin_rescored,  // PIN (MSBooster-edited or raw)
        ch_search_pepxml,        // pepXML (CrystalC-filtered or raw)
        ch_mzml_downstream,
        ch_fasta,
        ch_tool_configs,
        ch_peptideprophet_config,
        ch_data_type
    )

    //
    // STEP 6: ProteinProphet + Philosopher Filter (aggregate -> per-sample)
    //
    def aggregate_meta = [id: 'fragpipe']

    FRAGPIPE_INFERENCE(
        FRAGPIPE_VALIDATE.out.pepxml,
        ch_fasta,
        ch_tool_configs,
        aggregate_meta
    )

    //
    // STEP 7+: Post-inference analysis (parallel branches)
    // All post-inference steps run in parallel. FragPipe orders them sequentially
    // for single-machine OOM prevention, but in Nextflow the executor handles
    // resource scheduling. No true data dependencies between these steps
    // (verified against CmdIonquant.java, CmdPtmshepherd.java, CmdMBGMatch.java,
    // CmdSpecLibGen.java, CmdOPair.java -- none read each other's output).
    //

    //
    // STEP 7a: PTMShepherd (optional, aggregate)
    //
    FRAGPIPE_PTM(
        FRAGPIPE_INFERENCE.out.results_dir,
        FRAGPIPE_INFERENCE.out.protxml,
        ch_mzml_downstream,
        ch_tool_configs,
        aggregate_meta,
        ch_fasta
    )

    //
    // STEP 7b: Quantification - IonQuant / TMTIntegrator / FreeQuant (aggregate)
    //
    FRAGPIPE_QUANT(
        FRAGPIPE_INFERENCE.out.results_dir,
        ch_mzml_downstream,
        ch_annotation,
        ch_fasta,
        ch_tool_configs,
        ch_tmtintegrator_config,
        aggregate_meta
    )

    //
    // STEP 7c: Spectral library generation (optional, aggregate)
    //
    FRAGPIPE_SPECLIB(
        FRAGPIPE_INFERENCE.out.psms,
        FRAGPIPE_INFERENCE.out.peptides,
        ch_mzml_downstream,
        ch_fasta,
        ch_tool_configs,
        aggregate_meta
    )

    //
    // STEP 7d: Glycoproteomics analysis (optional, per-sample or aggregate)
    //
    // Glycan databases are only needed when MBG/OPair are enabled.
    // Use Channel.empty() so FRAGPIPE_GLYCO's maps produce empty channels
    // and glyco modules never get invoked when not configured.
    ch_glycan_dbs = Channel.empty()

    FRAGPIPE_GLYCO(
        FRAGPIPE_INFERENCE.out.psms,
        ch_mzml_downstream,
        ch_manifest,
        ch_glycan_dbs,
        ch_tool_configs
    )

    //
    // STEP 7e: Export - Skyline, SAINT, FPOP, Metaproteomics (optional)
    //
    // Taxon and SAINT files are only needed when their tools are enabled.
    // Use Channel.empty() so FRAGPIPE_EXPORT's maps produce empty channels
    // and export modules never get invoked when not configured.
    ch_taxon_files = Channel.empty()
    ch_saint_files = Channel.empty()

    // Speclib may be empty if FRAGPIPE_SPECLIB was disabled.
    // FRAGPIPE_EXPORT handles empty speclib with its own .ifEmpty(file('NO_SPECLIB')) sentinel.
    FRAGPIPE_EXPORT(
        FRAGPIPE_INFERENCE.out.results_dir,
        ch_mzml_downstream,
        FRAGPIPE_SPECLIB.out.library_speclib,
        ch_fasta,
        ch_taxon_files,
        ch_saint_files,
        ch_tool_configs,
        aggregate_meta
    )

    emit:
    // Converted spectra
    mzml                  = ch_mzml_downstream                                 // channel: [ val(meta), path(mzml) ]

    // Search outputs
    pepxml                = ch_search_pepxml                                   // channel: [ val(meta), path(pepxml) ]
    pin                   = ch_search_pin                                      // channel: [ val(meta), path(pin) ]

    // Validation outputs
    validated_pepxml      = FRAGPIPE_VALIDATE.out.pepxml                       // channel: [ val(meta), path(pepxml) ]

    // Inference outputs
    protxml               = FRAGPIPE_INFERENCE.out.protxml                     // channel: [ val(meta), path(combined.prot.xml) ]
    psms                  = FRAGPIPE_INFERENCE.out.psms                        // channel: [ val(meta), path(psm.tsv) ]
    peptides              = FRAGPIPE_INFERENCE.out.peptides                    // channel: [ val(meta), path(peptide.tsv) ]
    proteins              = FRAGPIPE_INFERENCE.out.proteins                    // channel: [ val(meta), path(protein.tsv) ]
    results_dir           = FRAGPIPE_INFERENCE.out.results_dir                 // channel: [ val(meta), path(results_dir) ]

    // Quantification outputs
    combined_protein      = FRAGPIPE_QUANT.out.combined_protein                // channel: [ val(meta), path(combined_protein.tsv) ]
    combined_peptide      = FRAGPIPE_QUANT.out.combined_peptide                // channel: [ val(meta), path(combined_peptide.tsv) ]
    combined_ion          = FRAGPIPE_QUANT.out.combined_ion                    // channel: [ val(meta), path(combined_ion.tsv) ]
    ionquant_dir          = FRAGPIPE_QUANT.out.ionquant_dir                    // channel: [ val(meta), path(results_dir) ]
    tmt_abundance         = FRAGPIPE_QUANT.out.tmt_abundance                   // channel: [ val(meta), path(abundance_*.tsv) ]
    tmt_dir               = FRAGPIPE_QUANT.out.tmt_dir                         // channel: [ val(meta), path(tmt-report) ]

    // PTM outputs
    ptm_results           = FRAGPIPE_PTM.out.results_dir                       // channel: [ val(meta), path(results_dir) ]

    // Spectral library outputs
    speclib               = FRAGPIPE_SPECLIB.out.library_speclib               // channel: [ val(meta), path(library.speclib) ]

    // Glyco outputs
    mbg_results           = FRAGPIPE_GLYCO.out.mbg_results                     // channel: [ val(meta), path(*_mbg*.tsv) ]
    opair_results         = FRAGPIPE_GLYCO.out.opair_results                   // channel: [ val(meta), path(*_opair_results.tsv) ]

    // Export outputs
    skyline_doc           = FRAGPIPE_EXPORT.out.skyline_doc                    // channel: [ val(meta), path(fragpipe.sky) ]
    taxonomy_results      = FRAGPIPE_EXPORT.out.taxonomy_results               // channel: [ val(meta), path(results_dir) ]

    // Versions collected via topic channels
    versions              = Channel.empty()
}
