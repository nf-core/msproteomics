/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_TOOL_CONFIGS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The tool configs the FragPipe subworkflows run from, by one of two planners (--fragpipe_planner):

    - parser (default): PARSEFRAGPIPEWORKFLOW reads the .workflow file.
    - plan: FRAGPIPEPLAN runs FragPipe headless --write-plan on the run's workflow,
      LC-MS files and database (a FragPipe image with --write-plan and --msfragger-per-file, passed with --fragpipe_container), and
      PLANTOTOOLCONFIGS turns FragPipe's work dir (run log, config files, run files) into the same outputs.

    Both emit the same channels, so the planner and the parser can be compared on every run.
    The plan needs the searched FASTA before the run, so --add_philosopher_decoy /
    --add_philosopher_contaminants (FragPipe's database update) work with the parser only
    (checkFragpipePlanner stops such a run in PIPELINE_INITIALISATION).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARSEFRAGPIPEWORKFLOW } from '../../../modules/local/parsefragpipeworkflow/main'
include { FRAGPIPEPLAN           } from '../../../modules/local/fragpipeplan/main'
include { PLANTOTOOLCONFIGS    } from '../../../modules/local/plantotoolconfigs/main'
include { parseFragpipeWorkflowTemplates } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

workflow FRAGPIPE_TOOL_CONFIGS {
    take:
    ch_file            // channel: [ val(sample_meta), path(spectra) ] - LC-MS files (meta: data_type, experiment, bioreplicate, group)
    ch_workflow_file   // channel: [ val(meta), path(workflow_file) ] - FragPipe workflow file
    ch_decoy_db        // channel: [ val(meta), path(fasta) ] - FASTA searched as given (empty with FragPipe's database update)
    ch_annotation      // channel: path(tmt_annotation.tsv) - value channel, [] when TMT-Integrator does not run
    val_planner        // string: 'parser' (PARSEFRAGPIPEWORKFLOW) or 'plan' (FRAGPIPEPLAN + PLANTOTOOLCONFIGS)

    main:
    if (val_planner == 'plan') {
        // FragPipe's manifest: one row per LC-MS file (headlessLayout), FragPipe's group from the samplesheet
        ch_rows = ch_file
            .map { meta, spectra -> [spectra.name, meta.experiment, meta.bioreplicate, meta.data_type, meta.group] }
            .toList()
        // TMT annotation rows (plex, channel, sample_name), as workflows/fragpipe.nf gives FRAGPIPE_HEADLESS_WF
        ch_annotation_rows = ch_annotation
            .flatMap { tsv -> tsv ? [tsv] : [] }
            .splitCsv(sep: '\t', header: true)
            .map { row -> [row.plex, row.channel, row.sample_name] }
            .toList()
        // original path of every staged input, for the portable plan, workflow and manifest (s3:// kept)
        ch_sources = ch_file
            .map { _meta, spectra -> [spectra.name, spectra.toAbsolutePath().normalize().toUriString()] }
            .toList()
            .map { pairs -> pairs.collectEntries { pair -> [pair[0], pair[1]] } }
        FRAGPIPEPLAN(
            ch_workflow_file
                .combine(ch_file.map { _meta, spectra -> spectra }.toList().map { files -> [files] })
                .combine(ch_decoy_db.map { _meta, fasta -> fasta })
                .combine(ch_rows.map { rows -> [rows] })
                .combine(ch_annotation_rows.map { rows -> [rows] })
                .combine(ch_sources.map { sources -> [sources] })
                .map { meta, workflow_file, files, fasta, rows, annotation_rows, sources ->
                    [meta, files, fasta, workflow_file.text, rows, annotation_rows, sources, fasta.toAbsolutePath().normalize().toUriString()]
                }
        )
        PLANTOTOOLCONFIGS(FRAGPIPEPLAN.out.workdir)
        ch_tool_configs_json       = PLANTOTOOLCONFIGS.out.tool_configs_json
        ch_msfragger_config        = PLANTOTOOLCONFIGS.out.msfragger_config
        ch_msfragger_dia_config    = PLANTOTOOLCONFIGS.out.msfragger_dia_config
        ch_msfragger_gpfdia_config = PLANTOTOOLCONFIGS.out.msfragger_gpfdia_config
        ch_msfragger_dda_plus      = PLANTOTOOLCONFIGS.out.msfragger_dda_plus_config
        ch_msbooster_config        = PLANTOTOOLCONFIGS.out.msbooster_config
        ch_portable_log            = FRAGPIPEPLAN.out.portable_log
        ch_portable_workflow       = FRAGPIPEPLAN.out.portable_workflow
        ch_portable_manifest       = FRAGPIPEPLAN.out.portable_manifest
    } else {
        // FragPipe adjusts the tools to the data types of all manifest files, DIA-Quant included
        // (TabWorkflow.adjustToolsBasedOnDataTypes), so the parser gets them with the workflow.
        ch_data_types = ch_file
            .map { meta, _spectra -> meta.data_type }
            .unique()
            .toSortedList()
        PARSEFRAGPIPEWORKFLOW(ch_workflow_file.combine(ch_data_types.map { data_types -> [data_types] }), parseFragpipeWorkflowTemplates(projectDir))
        ch_tool_configs_json       = PARSEFRAGPIPEWORKFLOW.out.tool_configs_json
        ch_msfragger_config        = PARSEFRAGPIPEWORKFLOW.out.msfragger_config
        ch_msfragger_dia_config    = PARSEFRAGPIPEWORKFLOW.out.msfragger_dia_config
        ch_msfragger_gpfdia_config = PARSEFRAGPIPEWORKFLOW.out.msfragger_gpfdia_config
        ch_msfragger_dda_plus      = PARSEFRAGPIPEWORKFLOW.out.msfragger_dda_plus_config
        ch_msbooster_config        = PARSEFRAGPIPEWORKFLOW.out.msbooster_config
        ch_portable_log            = channel.empty()
        ch_portable_workflow       = channel.empty()
        ch_portable_manifest       = channel.empty()
    }

    emit:
    tool_configs_json         = ch_tool_configs_json        // channel: [ val(meta), path(tool_configs.json) ]
    msfragger_config          = ch_msfragger_config         // channel: [ val(meta), path(fragger.params | fragger_dda.params) ]
    msfragger_dia_config      = ch_msfragger_dia_config     // channel: [ val(meta), path(fragger_dia.params) ]
    msfragger_gpfdia_config   = ch_msfragger_gpfdia_config  // channel: [ val(meta), path(fragger_gpfdia.params) ]
    msfragger_dda_plus_config = ch_msfragger_dda_plus       // channel: [ val(meta), path(fragger_dda_plus.params) ]
    msbooster_config          = ch_msbooster_config         // channel: [ val(meta), path(msbooster.config) ]
    fragpipe_log              = ch_portable_log             // channel: [ val(meta), path(log_*.txt) ] - plan mode: FragPipe's run log, licence redacted, original input paths, @WORKDIR@, @TOOLS_DIR@
    fragpipe_workflow         = ch_portable_workflow        // channel: [ val(meta), path(fragpipe.workflow) ] - plan mode: FragPipe's saved workflow, the same way
    fragpipe_manifest         = ch_portable_manifest        // channel: [ val(meta), path(fragpipe-files.fp-manifest) ] - plan mode: FragPipe's saved manifest
}
