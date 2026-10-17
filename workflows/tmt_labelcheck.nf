/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TMT LABEL CHECK WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A label check workflow is an ordinary FragPipe workflow (assets/TMT-labelcheck-{229,304}.workflow:
    TMT as a variable modification). MSPROTEOMICS_FRAGPIPE runs it, in pipeline or headless mode, with
    every FragPipe check; TMTLABELCHECKANALYZE then computes the labeling efficiency from each group's
    Report psm.tsv and, when IonQuant runs, from combined_modified_peptide.tsv.

    Labeling efficiency:
        Total_Sites   = peptide.count('K') + 1  (lysines + N-terminus)
        Labeled_Sites = count(TMT modifications)
        Efficiency    = sum(Labeled_Sites) / sum(Total_Sites) * 100
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MSPROTEOMICS_FRAGPIPE } from './fragpipe'
include { TMTLABELCHECKANALYZE                                } from '../modules/local/tmtlabelcheckanalyze/main'
include { TMTLABELCHECKANALYZE as TMTLABELCHECKANALYZE_IONQUANT } from '../modules/local/tmtlabelcheckanalyze/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MSPROTEOMICS_TMT_LABELCHECK {

    main:
    // --tmt_type names the label TMTLABELCHECKANALYZE counts, and selects the shipped label check workflow
    // (MSPROTEOMICS_FRAGPIPE; the label check takes no --fragpipe_workflow)
    if (!params.tmt_type) {
        error("--tmt_mode labelcheck requires --tmt_type (e.g. --tmt_type TMT16): the label the label check counts")
    }

    MSPROTEOMICS_FRAGPIPE()

    // Per-group labeling efficiency from each group's Report psm.tsv
    TMTLABELCHECKANALYZE(
        MSPROTEOMICS_FRAGPIPE.out.group_psms.map { psms -> [[id: 'tmt_labelcheck'], psms] },
        params.tmt_type,
        'psm'
    )

    // Combined labeling efficiency from IonQuant's combined_modified_peptide.tsv
    TMTLABELCHECKANALYZE_IONQUANT(
        MSPROTEOMICS_FRAGPIPE.out.combined_modified_peptide,
        params.tmt_type,
        'ionquant'
    )

    emit:
    multiqc_files = MSPROTEOMICS_FRAGPIPE.out.multiqc_files
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
