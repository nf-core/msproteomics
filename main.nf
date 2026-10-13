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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MSPROTEOMICS_DIANN          } from './workflows/diann'
include { MSPROTEOMICS_TMT_LABELCHECK } from './workflows/tmt_labelcheck'
include { MSPROTEOMICS_FRAGPIPE       } from './workflows/fragpipe'
include { MULTIQC                     } from './modules/nf-core/multiqc/main'
include { PIPELINE_INITIALISATION     } from './subworkflows/local/utils_nfcore_msproteomics_pipeline'
include { PIPELINE_COMPLETION         } from './subworkflows/local/utils_nfcore_msproteomics_pipeline'
include { methodsDescriptionText      } from './subworkflows/local/utils_nfcore_msproteomics_pipeline'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from './subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Route on --mode / --tmt_mode, then collate versions and run MultiQC.
// Each mode workflow parses the samplesheet itself (samplesheetToList with
// assets/schema_input.json, see parseSamplesheet), because FragPipe's group rule
// depends on the mode's .workflow file.
//
workflow NFCORE_MSPROTEOMICS {

    take:
    outdir // string: output directory

    main:
    if (params.mode == 'diann') {
        MSPROTEOMICS_DIANN()
        ch_versions      = MSPROTEOMICS_DIANN.out.versions
        ch_multiqc_files = MSPROTEOMICS_DIANN.out.multiqc_files
    } else if (params.mode == 'fragpipe' && params.tmt_mode == 'labelcheck') {
        MSPROTEOMICS_TMT_LABELCHECK()
        ch_versions      = MSPROTEOMICS_TMT_LABELCHECK.out.versions
        ch_multiqc_files = MSPROTEOMICS_TMT_LABELCHECK.out.multiqc_files
    } else if (params.mode == 'fragpipe') {
        MSPROTEOMICS_FRAGPIPE()
        ch_versions      = MSPROTEOMICS_FRAGPIPE.out.versions
        ch_multiqc_files = MSPROTEOMICS_FRAGPIPE.out.multiqc_files
    } else {
        error("--mode must be 'diann' or 'fragpipe', got '${params.mode}'")
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_' + 'msproteomics_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def multiqc_configs = [file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)] +
        (params.multiqc_config ? [file(params.multiqc_config, checkIfExists: true)] : [])
    def multiqc_logo = params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : []
    def methods_description_file = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)

    def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
        .collectFile(name: 'workflow_summary_mqc.yaml')
    def ch_methods_description = channel.value(methodsDescriptionText(methods_description_file))
        .collectFile(name: 'methods_description_mqc.yaml', sort: true)

    MULTIQC(
        ch_multiqc_files
            .mix(ch_workflow_summary, ch_collated_versions, ch_methods_description)
            .collect()
            .map { files -> [[id: 'msproteomics'], files, multiqc_configs, multiqc_logo, [], []] }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> report } // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Version, help, parameter validation (nextflow_schema.json) and config checks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_MSPROTEOMICS(params.outdir)

    //
    // SUBWORKFLOW: Completion email and summary
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_MSPROTEOMICS.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
