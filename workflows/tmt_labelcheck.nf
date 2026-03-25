/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TMT LABEL CHECK WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_msproteomics_pipeline'

// SDRF generation (bookkeeping only)
include { GENERATE_SDRF_FROM_SAMPLESHEET } from '../modules/local/generate_sdrf_from_samplesheet/main'

// TMT Label Check subworkflow
include { TMT_LABELCHECK } from '../subworkflows/local/tmt_labelcheck/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MSPROTEOMICS_TMT_LABELCHECK {

    take:
    ch_input    // channel: [ meta, spectra_file ]

    main:

    GENERATE_SDRF_FROM_SAMPLESHEET(
        Channel.fromPath(params.input, checkIfExists: true),
        WorkflowUtils.generateSdrfBookkeeping(params)
    )

    def organism = WorkflowUtils.convertOrganismToStandardName(params.organism ?: 'Homo sapiens')
    def db_channels = WorkflowUtils.resolveFragPipeDatabaseChannels(organism, params)

    // Auto-select TMT labelcheck workflow file based on tmt_type if not explicitly provided
    def tmt_workflow = params.fragpipe_workflow ?: WorkflowUtils.getTmtLabelcheckWorkflow(params.tmt_type, projectDir)
    ch_workflow_file = Channel.of(
        [[id: 'labelcheck_workflow'], file(tmt_workflow, checkIfExists: true)]
    )

    TMT_LABELCHECK(
        ch_input,
        db_channels.ch_fasta,
        ch_workflow_file,
        params.tmt_type,
        db_channels.ch_prebuilt_db
    )

    emit:
    multiqc_report = Channel.empty()
    versions = Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
