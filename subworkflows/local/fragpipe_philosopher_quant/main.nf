/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_PHILOSOPHER_QUANT SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The per-group tables after Philosopher Filter, in FragPipe 24.0's order (FragpipeRun.java:2463-2471,
    :2490): FREEQUANT (label-free, or the TMT FreeQuant when the TMT-Integrator extraction tool is
    Philosopher) -> PHILOSOPHERLABELQUANT (TMT) -> PHILOSOPHERREPORT on the workspace Filter left
    (groups with meta.report_after_quant), then ExportMatchedFragments on every group.

    Its output is what every later tool reads: the glyco chain (PairScans, O-Pair, PTM-Shepherd, MBG
    have IN-edges from PhilosopherReport), IonQuant, TMT-Integrator, SpecLibGen and the exports.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FREEQUANT                } from '../../../modules/local/freequant/main'
include { PHILOSOPHERLABELQUANT } from '../../../modules/local/philosopherlabelquant/main'
include { WRITETEXTFILE as WRITE_LABELQUANT_ANNOTATION } from '../../../modules/local/writetextfile/main'
include { WRITETEXTFILE as WRITE_LABELQUANT_EXPERIMENT_ANNOTATION } from '../../../modules/local/writetextfile/main'
include { labelquantAnnotationText } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'
include { PHILOSOPHERREPORT } from '../../../modules/local/philosopherreport/main'
include { shouldRunTool; getToolArgs; getToolField; getToolReportArgs; groupResultsDir; isobaricExperimentAnnotationLines } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'
include { FRAGPIPE_MATCHED_FRAGMENTS } from '../fragpipe_matched_fragments/main'

workflow FRAGPIPE_PHILOSOPHER_QUANT {
    take:
    ch_results_dir          // channel: [ val(meta), path(results_dir) ] - per-group PHILOSOPHERFILTER output (meta.group, meta.report_after_quant)
    ch_mzml                 // channel: [ val(meta), path(mzml) ] - searched spectra (meta.group = FragPipe group)
    ch_lcms_files           // channel: [ val(meta), path(spectra) ] - every samplesheet LC-MS file as given (FragPipe's file list)
    ch_annotation           // channel: path(annotation_file) - multi-plex TMT annotation, [] when TMT-Integrator does not run
    ch_tool_configs         // channel: val(tool_configs_map) - parsed JSON config
    ch_emf_inputs           // channel: [ path(manifest), [spectra], [umpire_spectra] ] - FRAGPIPE_MATCHED_FRAGMENTS inputs
    aggregate_meta          // val(meta) for aggregate outputs

    main:

    // FragPipe's LC-MS file list (wrapped so combine keeps it one item) and the TMT annotation
    // lines (header first; [] without TMT)
    ch_lcms_list = ch_lcms_files.toList().map { files -> [files] }
    // ch_annotation is a value channel: the annotation file is read once
    ch_annotation_lines = ch_annotation.map { annot -> annot ? [annot.readLines()] : [[]] }

    //
    // STEP 1: FreeQuant per group (label-free, FragpipeRun.java:1812-1820, or the TMT FreeQuant
    // before LabelQuant, FragpipeRun.java:2049-2055), in the group dir Filter left. FreeQuant and
    // LabelQuant take mzML only (CmdFreequant.java:44, CmdLabelquant.java:51); every LC-MS file
    // of a group is passed in one --dir.
    //
    ch_group_lcms = ch_mzml
        .map { meta, mzml -> [meta.group, mzml] }
        // by file name: remote (http/S3) and task-dir paths of one group do not compare with each other
        .groupTuple(sort: { a, b -> a.name <=> b.name })

    ch_philosopher_quant = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'freequant') || shouldRunTool(configs, 'labelquant') }
        .combine(ch_lcms_list)
        .map { configs, lcms_files ->
            def unsupported = lcms_files.findAll { item -> item[0].data_type != 'DIA-Quant' }
                .collect { item -> item[1].name.toLowerCase().tokenize('.').last() }
                .unique().findAll { ext -> ext != 'mzml' }
            if (unsupported) {
                error(shouldRunTool(configs, 'labelquant') ?
                    "Philosopher as the TMT intensity extraction tool doesn't support '.${unsupported.join(', ')}' files. Please replace with mzML format or set IonQuant as the TMT intensity extraction tool." :
                    "FreeQuant doesn't support '.${unsupported.join(', ')}' files. Please switch to IonQuant in the Quant (MS1) tab.")
            }
            configs
        }

    // Both sides only when FreeQuant / LabelQuant run, so every group has its Filter dir
    ch_phi_group_lcms = ch_group_lcms
        .combine(ch_philosopher_quant)
        .map { group, lcms_files, _configs -> [group, lcms_files] }

    ch_for_freequant = ch_results_dir
        .combine(ch_philosopher_quant)
        .map { meta, dir, configs -> [meta.group, meta, dir, configs] }
        .join(ch_phi_group_lcms, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .map { _group, meta, dir, configs, lcms_files -> [meta, dir, lcms_files, getToolArgs(configs, 'freequant')] }

    FREEQUANT(ch_for_freequant)

    //
    // STEP 2: Philosopher LabelQuant per group (CmdLabelquant.java:62-170). Each group's <group>_annotation.txt
    // comes from the multi-plex TMT annotation (labelquantAnnotationText); the isobaric experiment_annotation.tsv
    // FragPipe writes once in the work dir is written with the first group (CmdLabelquant.java:163-170).
    //
    ch_labelquant_groups = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'labelquant') }
        .combine(
            FREEQUANT.out.results_dir
                .map { meta, dir -> [meta.group, meta, dir] }
                .join(ch_phi_group_lcms, by: 0, failOnMismatch: true, failOnDuplicate: true)
        )
        .combine(ch_annotation_lines)
        .combine(ch_results_dir.map { meta, _dir -> meta.group }.collect(sort: true).map { groups -> [groups] })

    WRITE_LABELQUANT_ANNOTATION(
        ch_labelquant_groups.map { _configs, group, meta, _dir, _lcms_files, annotation_lines, _groups ->
            [meta, "${group}_annotation.txt".toString(), labelquantAnnotationText(annotation_lines, group)]
        }
    )

    WRITE_LABELQUANT_EXPERIMENT_ANNOTATION(
        ch_labelquant_groups
            .filter { _configs, group, _meta, _dir, _lcms_files, _annotation_lines, groups -> group == groups[0] }
            .map { _configs, _group, meta, _dir, _lcms_files, annotation_lines, _groups ->
                [meta, 'experiment_annotation.tsv', isobaricExperimentAnnotationLines(annotation_lines).collect { line -> "${line}\n" }.join('')]
            }
            .filter { _meta, _name, text -> text }
    )

    ch_for_labelquant = ch_labelquant_groups
        .map { configs, _group, meta, dir, lcms_files, _annotation_lines, _groups ->
            [meta, dir, lcms_files, getToolArgs(configs, 'labelquant'), getToolField(configs, 'labelquant', 'brand')]
        }
        .join(WRITE_LABELQUANT_ANNOTATION.out.file, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, dir, lcms_files, config_cli, brand, annotation -> [meta, dir, lcms_files, annotation, config_cli, brand] }

    PHILOSOPHERLABELQUANT(ch_for_labelquant)

    //
    // STEP 3: Philosopher Report per group after FreeQuant / LabelQuant (FragpipeRun.java:2471),
    // then the workspace clean
    //
    ch_for_report = PHILOSOPHERLABELQUANT.out.results_dir
        .mix(
            FREEQUANT.out.results_dir
                .combine(ch_tool_configs)
                .filter { _meta, _dir, configs -> !shouldRunTool(configs, 'labelquant') }
                .map { meta, dir, _configs -> [meta, dir] }
        )
        .combine(ch_tool_configs)
        .map { meta, dir, configs -> [meta, dir, getToolReportArgs(configs, 'filter')] }

    PHILOSOPHERREPORT(ch_for_report)

    ch_report_results_dir = PHILOSOPHERREPORT.out.results_files.map { meta, files -> [meta, groupResultsDir(files)] }

    // Final per-group tables every later reader sees: the Report output after FreeQuant/LabelQuant,
    // else the Filter output (Filter leaves only placeholder tables when Report follows, ext.skip_report),
    // with matched fragments exported right after the Report (FragpipeRun.java:2490)
    // FragPipe exports matched fragments right after the final PhilosopherReport and before every
    // tool that reads psm.tsv (FragpipeRun.java:2411-2417, :2490); it rewrites psm.tsv in place.
    FRAGPIPE_MATCHED_FRAGMENTS(
        ch_results_dir
            .filter { meta, _dir -> !meta.report_after_quant }
            .mix(ch_report_results_dir),
        ch_emf_inputs,
        ch_tool_configs,
        aggregate_meta
    )
    ch_final_results = FRAGPIPE_MATCHED_FRAGMENTS.out.results_dir

    emit:
    final_results_dir      = ch_final_results                         // channel: [ val(meta), path(results_dir) ] - per-group tables after the final report and the matched-fragment export
    experiment_annotation  = WRITE_LABELQUANT_EXPERIMENT_ANNOTATION.out.file  // channel: [ val(meta), path(experiment_annotation.tsv) ]
}
