/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_VALIDATE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PSM validation and scoring refinement.

    Modules:
    - PERCOLATOR (machine learning PSM rescoring)
    - PEPTIDEPROPHET (statistical PSM validation, optional alternative to Percolator)
    - PTMPROPHET (PTM site localization, optional)

    Execution: Per-sample parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PERCOLATOR     } from '../../../modules/local/percolator/main'
include { PEPTIDEPROPHET } from '../../../modules/local/philosopher/peptideprophet/main'
include { PTMPROPHET     } from '../../../modules/local/ptmprophet/main'
include { shouldRunTool; getToolArgs } from '../fragpipe_utils'

workflow FRAGPIPE_VALIDATE {
    take:
    ch_pin                    // channel: [ val(meta), path(pin) ] - PIN files from MSFragger/MSBooster
    ch_pepxml                 // channel: [ val(meta), path(pepxml) ] - pepXML from MSFragger
    ch_mzml                   // channel: [ val(meta), path(mzml) ] - mzML files
    ch_fasta                  // channel: path(fasta) - value channel for database
    ch_tool_configs           // channel: val(tool_configs_map) - parsed JSON config
    ch_peptideprophet_config  // channel: path(config_file) - value channel, [] when not used
    data_type                 // val: 'DDA' or 'DIA'

    main:

    //
    // STEP 1: Percolator PSM rescoring (per-sample)
    // Uses PIN file from MSFragger/MSBooster and outputs validated pepXML
    //
    ch_pin_keyed = ch_pin
        .map { meta, pin -> [meta.id, meta, pin] }

    ch_pepxml_keyed = ch_pepxml
        .map { meta, pepxml -> [meta.id, meta, pepxml] }

    ch_mzml_keyed = ch_mzml
        .map { meta, mzml -> [meta.id, mzml] }

    ch_for_percolator = ch_pin_keyed
        .join(ch_pepxml_keyed, by: 0, failOnMismatch: false, failOnDuplicate: false)
        .map { key, meta_pin, pin, meta_pepxml, pepxml -> [key, meta_pin, pin, pepxml] }
        .join(ch_mzml_keyed, by: 0, failOnMismatch: false, failOnDuplicate: false)
        .combine(ch_tool_configs)
        .filter { key, meta, pin, pepxml, mzml, configs -> shouldRunTool(configs, 'percolator') }
        .map { key, meta, pin, pepxml, mzml, configs ->
            def args = getToolArgs(configs, 'percolator')
            [meta, pin, pepxml, mzml, args]
        }

    PERCOLATOR(ch_for_percolator, data_type)


    //
    // STEP 2: PeptideProphet (optional, alternative to Percolator)
    // Used when Percolator is not run
    //
    ch_pepxml_for_peptideprophet = ch_pepxml_keyed
        .join(ch_mzml_keyed, by: 0, failOnMismatch: false, failOnDuplicate: false)
        .combine(ch_tool_configs)
        .filter { key, meta, pepxml, mzml, configs ->
            shouldRunTool(configs, 'peptideprophet') && !shouldRunTool(configs, 'percolator')
        }
        .map { key, meta, pepxml, mzml, configs ->
            [meta, pepxml, mzml]
        }
        .combine(ch_peptideprophet_config)
        .map { meta, pepxml, mzml, config -> [meta, pepxml, mzml, config] }

    PEPTIDEPROPHET(ch_pepxml_for_peptideprophet, ch_fasta)


    //
    // STEP 3: PTMProphet (optional, for PTM site localization)
    // Runs on pepXML output from Percolator or PeptideProphet
    //
    ch_pepxml_from_percolator = PERCOLATOR.out.pepxml
        .map { meta, pepxml -> [meta.id, meta, pepxml] }

    ch_pepxml_from_peptideprophet = PEPTIDEPROPHET.out.pepxml
        .map { meta, pepxml -> [meta.id, meta, pepxml] }

    ch_validated_pepxml = ch_pepxml_from_percolator
        .mix(ch_pepxml_from_peptideprophet)

    ch_for_ptmprophet = ch_validated_pepxml
        .combine(ch_tool_configs)
        .filter { key, meta, pepxml, configs -> shouldRunTool(configs, 'ptmprophet') }
        .map { key, meta, pepxml, configs ->
            def args = getToolArgs(configs, 'ptmprophet')
            [meta, pepxml, args]
        }

    PTMPROPHET(ch_for_ptmprophet)


    //
    // Collect outputs
    // pepXML: Use PTMProphet output if available, else validated pepXML, else passthrough.
    // Routing uses actual output presence (left-join pattern) rather than config flags alone,
    // so if PTMProphet is enabled but produces no output for a sample, that sample's
    // validated pepXML is still emitted.
    //
    ch_ptmprophet_keyed = PTMPROPHET.out.mod_pepxml
        .map { meta, pepxml -> [meta.id, pepxml] }

    // Left-join validated pepXML with PTMProphet output: prefer PTMProphet if available
    ch_pepxml_final = ch_validated_pepxml
        .join(ch_ptmprophet_keyed, by: 0, remainder: true)
        .map { key, meta, validated_pepxml, ptm_pepxml ->
            def pepxml = ptm_pepxml ?: validated_pepxml
            [meta, pepxml]
        }

    // Also output unvalidated pepXML for workflows that don't use Percolator/PeptideProphet
    ch_pepxml_passthrough = ch_pepxml_keyed
        .combine(ch_tool_configs)
        .filter { key, meta, pepxml, configs ->
            !shouldRunTool(configs, 'percolator') && !shouldRunTool(configs, 'peptideprophet')
        }
        .map { key, meta, pepxml, configs -> [meta, pepxml] }

    ch_pepxml_all = ch_pepxml_final.mix(ch_pepxml_passthrough)

    emit:
    pepxml       = ch_pepxml_all                              // channel: [ val(meta), path(pepxml) ] - validated or passthrough
    target_psms  = PERCOLATOR.out.target_psms                 // channel: [ val(meta), path(target_psms.tsv) ]
    decoy_psms   = PERCOLATOR.out.decoy_psms                  // channel: [ val(meta), path(decoy_psms.tsv) ]
    mod_pepxml   = PTMPROPHET.out.mod_pepxml                  // channel: [ val(meta), path(mod.pep.xml) ] - PTM localized
    results_dir  = PERCOLATOR.out.results_dir                 // channel: [ val(meta), path(results_dir) ]
    versions     = PERCOLATOR.out.versions_percolator
        .mix(PERCOLATOR.out.versions_fragpipe)
        .mix(PEPTIDEPROPHET.out.versions_philosopher)
        .mix(PEPTIDEPROPHET.out.versions_fragpipe)
        .mix(PTMPROPHET.out.versions_ptmprophet)
        .mix(PTMPROPHET.out.versions_fragpipe)
}
