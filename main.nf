#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/msproteomics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/msproteomics
    Website: https://nf-co.re/msproteomics
    Slack  : https://nfcore.slack.com/channels/msproteomics
----------------------------------------------------------------------------------------
*/

include { MSPROTEOMICS_DIANN          } from './workflows/diann'
include { MSPROTEOMICS_TMT_LABELCHECK } from './workflows/tmt_labelcheck'
include { MSPROTEOMICS_FRAGPIPE       } from './workflows/fragpipe'
include { PIPELINE_INITIALISATION     } from './subworkflows/local/utils_nfcore_msproteomics_pipeline'
include { PIPELINE_COMPLETION         } from './subworkflows/local/utils_nfcore_msproteomics_pipeline'

workflow {
    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow based on --mode
    //
    if (params.mode == 'diann') {
        MSPROTEOMICS_DIANN(PIPELINE_INITIALISATION.out.samplesheet)
    } else if (params.mode == 'fragpipe') {
        if (params.tmt_mode == 'labelcheck') {
            MSPROTEOMICS_TMT_LABELCHECK(PIPELINE_INITIALISATION.out.samplesheet)
        } else {
            MSPROTEOMICS_FRAGPIPE(PIPELINE_INITIALISATION.out.samplesheet)
        }
    } else {
        error "Parameter 'mode' is required. Use --mode diann or --mode fragpipe"
    }

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        Channel.empty()
    )
}
