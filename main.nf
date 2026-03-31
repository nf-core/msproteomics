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

workflow {
    if (params.mode == 'diann') {
        MSPROTEOMICS_DIANN()
    } else if (params.mode == 'fragpipe') {
        if (params.tmt_mode == 'labelcheck') {
            MSPROTEOMICS_TMT_LABELCHECK()
        } else {
            MSPROTEOMICS_FRAGPIPE()
        }
    } else {
        error "Parameter 'mode' is required. Use --mode diann or --mode fragpipe"
    }
}
