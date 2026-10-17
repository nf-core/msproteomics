/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_QUANT SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Quantification after the glyco chain, in FragPipe 24.0's command order and data flow
    (FragpipeRun.java:2477-2481). The Philosopher path (FreeQuant, LabelQuant, final Report) runs
    before the glyco chain in FRAGPIPE_PHILOSOPHER_QUANT.

    - IONQUANT (label-free) -> IONQUANT_MS1 (TMT) -> IONQUANT_ISOBARIC (TMT): FragPipe runs
      cmdIonquant, cmdTmtIonquant and cmdTmtIonquantIsobaric one after another in the same
      work dir, so each pass reads the psm.tsv files the previous pass rewrote
    - TMTINTEGRATOR on one psm.tsv per group (IonQuant isobaric output, or the final tables of
      the Philosopher TMT path)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IONQUANT                      } from '../../../modules/local/ionquant/main'
include { IONQUANT as IONQUANT_MS1      } from '../../../modules/local/ionquant/main'
include { IONQUANT as IONQUANT_ISOBARIC } from '../../../modules/local/ionquant/main'
include { TMTINTEGRATOR            } from '../../../modules/local/tmtintegrator/main'
include { TMTINTEGRATOR as TMTINTEGRATOR_UNMOD } from '../../../modules/local/tmtintegrator/main'
include { WRITETEXTFILE as WRITE_MSSTATS_TMT_ANNOTATION } from '../../../modules/local/writetextfile/main'
include { shouldRunTool; getToolArgs; getToolField; getToolModmasses; fragpipeRootGroupDir; ionquantGroupDirs; lfqExperimentAnnotationLines; isobaricExperimentAnnotationLines; msstatsTmtAnnotationLines; ionquantInputDataTypeArgs; fpopCoadaptrMeta } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

workflow FRAGPIPE_QUANT {
    take:
    ch_results_dir          // channel: [ val(meta), path(results_dir) ] - final per-group tables (FRAGPIPE_PHILOSOPHER_QUANT, then the glyco chain; meta.group)
    ch_mzml                 // channel: [ val(meta), path(mzml) ] - searched spectra (meta.group = FragPipe group)
    ch_lcms_files           // channel: [ val(meta), path(spectra) ] - every samplesheet LC-MS file as given (FragPipe's file list)
    ch_annotation           // channel: path(annotation_file) - multi-plex TMT annotation, [] when TMT-Integrator does not run
    ch_tool_configs         // channel: val(tool_configs_map) - parsed JSON config
    ch_glyco_masses         // channel: path(glyco_masses_list.txt) - value channel, [] when AppendFile does not run (FragpipeRun.java:1822-1830)
    aggregate_meta          // val(meta) for aggregate outputs

    main:

    // FragPipe's LC-MS file list (wrapped so combine keeps it one item) and the TMT annotation
    // lines (header first; [] without TMT)
    ch_lcms_list = ch_lcms_files.toList().map { files -> [files] }
    // ch_annotation is a value channel: the annotation file is read once
    ch_annotation_lines = ch_annotation.map { annot -> annot ? [annot.readLines()] : [[]] }

    // Per-group results as [ [group dirs], root_group_dir ]: root_group_dir names the blank
    // group "" ('' if none), which IonQuant reads from its work-dir root like FragPipe.
    ch_all_results = ch_results_dir
        .collect(flat: false)
        .map { items ->
            def sorted = items.toSorted { item -> item[0].id }
            [sorted.collect { item -> item[1] }, fragpipeRootGroupDir(sorted.collect { item -> item[0] })]
        }

    ch_all_mzml = ch_mzml
        .map { _meta, mzml -> mzml }
        .collect(sort: true)
        .map { items -> [items] }  // Wrap in list to prevent combine from flattening

    //
    // STEP 4: label-free IonQuant (FragpipeRun.java:1859-1891). It writes the label-free
    // experiment_annotation.tsv unless a TMT pass configured after it in FragPipe replaces the
    // file (IonQuant isobaric pass, FragpipeRun.java:2005-2028; LabelQuant, CmdLabelquant.java:163-170).
    //
    ch_for_ionquant = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'ionquant') }
        .combine(ch_all_results)
        .combine(ch_all_mzml)
        .combine(ch_lcms_list)
        .map { configs, results, root_group_dir, mzml, lcms_files ->
            def tmt_annotation_later = shouldRunTool(configs, 'ionquant_isobaric') || shouldRunTool(configs, 'labelquant')
            [
                fpopCoadaptrMeta(aggregate_meta, configs, ['ionquant']), results, mzml,
                ionquantInputDataTypeArgs(getToolArgs(configs, 'ionquant'), lcms_files),
                getToolModmasses(configs, 'ionquant'), root_group_dir,
                tmt_annotation_later ? [] : lfqExperimentAnnotationLines(lcms_files)
            ]
        }

    IONQUANT(
        ch_for_ionquant.map { item -> [item[0], item[1], item[2]] },
        [],
        ch_for_ionquant.map { item -> item[3] },
        ch_for_ionquant.map { item -> item[4] },
        ch_for_ionquant.map { item -> item[5] },
        ch_for_ionquant.map { item -> item[6] },
        ch_glyco_masses
    )

    //
    // STEP 5: TMT IonQuant MS1 pass (FragpipeRun.java:1978-1995). It reads the label-free IonQuant
    // output when that ran before it, else the final per-group tables.
    //
    // root_group exists only with a blank group (an empty output directory has no object on S3):
    // the remainder join gives null for its absence.
    ch_ms1_results = IONQUANT.out.results_dir
        .join(IONQUANT.out.root_group, remainder: true, failOnDuplicate: true)
        .map { _meta, dir, root_dir -> ionquantGroupDirs(dir, root_dir) }
        .mix(
            ch_tool_configs
                .filter { configs -> !shouldRunTool(configs, 'ionquant') }
                .combine(ch_all_results)
                .map { _configs, results, root_group_dir -> [results, root_group_dir] }
        )

    ch_for_ionquant_ms1 = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'ionquant_ms1') }
        .combine(ch_ms1_results)
        .combine(ch_all_mzml)
        .combine(ch_lcms_list)
        .map { configs, results, root_group_dir, mzml, lcms_files ->
            [
                aggregate_meta, results, mzml,
                ionquantInputDataTypeArgs(getToolArgs(configs, 'ionquant_ms1'), lcms_files),
                getToolModmasses(configs, 'ionquant_ms1'), root_group_dir
            ]
        }

    IONQUANT_MS1(
        ch_for_ionquant_ms1.map { item -> [item[0], item[1], item[2]] },
        [],
        ch_for_ionquant_ms1.map { item -> item[3] },
        ch_for_ionquant_ms1.map { item -> item[4] },
        ch_for_ionquant_ms1.map { item -> item[5] },
        [],
        ch_glyco_masses
    )

    //
    // STEP 6: TMT IonQuant isobaric pass on the MS1-updated psm.tsv files (FragpipeRun.java:2005-2028),
    // writing the isobaric experiment_annotation.tsv
    //
    ch_for_ionquant_isobaric = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'ionquant_isobaric') }
        .combine(
            IONQUANT_MS1.out.results_dir
                .join(IONQUANT_MS1.out.root_group, remainder: true, failOnDuplicate: true)
                .map { _meta, dir, root_dir -> ionquantGroupDirs(dir, root_dir) }
        )
        .combine(ch_all_mzml)
        .combine(ch_lcms_list)
        .combine(ch_annotation_lines)
        .map { configs, results, root_group_dir, mzml, lcms_files, annotation_lines ->
            [
                fpopCoadaptrMeta(aggregate_meta, configs, ['isobaric']), results, mzml,
                ionquantInputDataTypeArgs(getToolArgs(configs, 'ionquant_isobaric'), lcms_files),
                getToolModmasses(configs, 'ionquant_isobaric'), root_group_dir,
                isobaricExperimentAnnotationLines(annotation_lines)
            ]
        }

    IONQUANT_ISOBARIC(
        ch_for_ionquant_isobaric.map { item -> [item[0], item[1], item[2]] },
        ch_annotation_lines.map { wrapped -> wrapped[0] },
        ch_for_ionquant_isobaric.map { item -> item[3] },
        ch_for_ionquant_isobaric.map { item -> item[4] },
        ch_for_ionquant_isobaric.map { item -> item[5] },
        ch_for_ionquant_isobaric.map { item -> item[6] },
        ch_glyco_masses
    )

    //
    // STEP 7: TMT-Integrator on one psm.tsv per group (CmdTmtIntegrator.java:288-296): the IonQuant
    // isobaric pass output, or the Philosopher report output. FragPipe skips TMT-Integrator when
    // any LC-MS file is DIA, DIA-Lib, DIA-Quant or GPF-DIA and accepts only mzML, raw and .d
    // (FragpipeRun.java:1896-1912).
    //
    ch_tmt_from_ionquant = IONQUANT_ISOBARIC.out.results_dir
        .join(IONQUANT_ISOBARIC.out.root_group, remainder: true, failOnDuplicate: true)
        .map { _meta, dir, root_dir -> ionquantGroupDirs(dir, root_dir) }

    ch_tmt_from_report = ch_results_dir
        .filter { meta, _dir -> meta.report_after_quant }
        .collect(flat: false)
        .map { items -> [items] }
        .combine(ch_tool_configs)
        .filter { _items, configs -> shouldRunTool(configs, 'labelquant') }
        .map { items, _configs ->
            def sorted = items.sort { item -> item[0].id }
            [sorted.collect { item -> item[1] }, fragpipeRootGroupDir(sorted.collect { item -> item[0] })]
        }

    ch_for_tmtintegrator = ch_tmt_from_ionquant
        .mix(ch_tmt_from_report)
        .combine(ch_tool_configs)
        .combine(ch_lcms_list)
        .combine(ch_annotation_lines)
        .filter { _dirs, _root, configs, lcms_files, _lines ->
            shouldRunTool(configs, 'tmtintegrator') &&
                !lcms_files.any { item -> item[0].data_type in ['DIA', 'DIA-Lib', 'DIA-Quant', 'GPF-DIA'] }
        }
        .map { dirs, root_group_dir, configs, lcms_files, annotation_lines ->
            if (lcms_files.any { item -> !(item[1].name.toLowerCase() ==~ /.+\.(mzml|raw|d)/) }) {
                error("TMT-Integrator only supports mzML, raw, and .d files.\nPlease remove other files from the input list.")
            }
            def msstats = getToolField(configs, 'tmtintegrator', 'philosopher_msstats', false) ?
                msstatsTmtAnnotationLines(annotation_lines, lcms_files) : []
            [aggregate_meta, dirs, root_group_dir, msstats]
        }

    // MSstatsTMT_annotation.csv (CmdTmtIntegrator.java:105-140), written once and copied into each run's work dir
    // The run's meta (aggregate_meta, one TMT-Integrator run per pipeline run) is the key that carries the CSV back into the run's tuple.
    WRITE_MSSTATS_TMT_ANNOTATION(
        ch_for_tmtintegrator
            .filter { _meta, _dirs, _root, msstats -> msstats }
            .map { meta, _dirs, _root, msstats -> [meta, 'MSstatsTMT_annotation.csv', msstats.collect { line -> "${line}\n" }.join('')] }
    )
    ch_tmti_runs = ch_for_tmtintegrator
        .map { meta, dirs, root_group_dir, _msstats -> [meta, dirs, root_group_dir] }
        .join(WRITE_MSSTATS_TMT_ANNOTATION.out.file, remainder: true, failOnDuplicate: true)
        .map { meta, dirs, root_group_dir, csv -> [meta, dirs, root_group_dir, csv ?: []] }

    TMTINTEGRATOR(
        ch_tmti_runs,
        // tmt-integrator-conf.yml text: FragPipe's plan tmtintegrator.args (tool_configs JSON, a channel value)
        ch_tool_configs.map { configs -> getToolArgs(configs, 'tmtintegrator') }
    )

    //
    //
    // TMT-Integrator second run for FPOP TMT (FragpipeRun.java:2195-2209): same PSM tables,
    // mod_tag none and min_site_prob -1, written to tmt-report-unmod (CmdTmtIntegrator.java:149-152).
    //
    TMTINTEGRATOR_UNMOD(
        ch_tmti_runs
            .combine(ch_tool_configs)
            .filter { row -> shouldRunTool(row[-1], 'tmtintegrator_unmod') }
            .map { meta, dirs, root_group_dir, csv, _configs -> [meta, dirs, root_group_dir, csv] },
        ch_tool_configs.map { configs -> getToolArgs(configs, 'tmtintegrator_unmod') }
    )

    // REPRINT tables IonQuant writes into the work dir (reprint.spc.tsv / reprint.int.tsv). FragPipe's
    // Downstream SAINTexpress reads the work dir's copy, i.e. the one of the last IonQuant run that wrote
    // it, in FragPipe order: LFQ IonQuant, or TMT MS1 then isobaric (FragpipeRun.java:2477-2479).
    //
    ch_reprint_runs = IONQUANT.out.reprint_spc.map { _meta, tsv -> ['spc', 0, tsv] }
        .mix(IONQUANT.out.reprint_int.map { _meta, tsv -> ['int', 0, tsv] })
        .mix(IONQUANT_MS1.out.reprint_spc.map { _meta, tsv -> ['spc', 1, tsv] })
        .mix(IONQUANT_MS1.out.reprint_int.map { _meta, tsv -> ['int', 1, tsv] })
        .mix(IONQUANT_ISOBARIC.out.reprint_spc.map { _meta, tsv -> ['spc', 2, tsv] })
        .mix(IONQUANT_ISOBARIC.out.reprint_int.map { _meta, tsv -> ['int', 2, tsv] })
    ch_reprint = ch_reprint_runs
        .toList()
        .flatMap { rows ->
            rows.groupBy { row -> row[0] }.collect { mode, mode_rows -> [mode, mode_rows.max { row -> row[1] }[2]] }
        }

    // FragPipe's work dir after the last IonQuant pass (label-free, TMT MS1, isobaric): each pass rewrites
    // the groups' psm.tsv and ion.tsv in place (FragpipeRun.java:2477-2479); empty when IonQuant does not run
    ch_quant_dir = IONQUANT.out.results_dir.map { _meta, dir -> [0, dir] }
        .mix(
            IONQUANT_MS1.out.results_dir.map { _meta, dir -> [1, dir] },
            IONQUANT_ISOBARIC.out.results_dir.map { _meta, dir -> [2, dir] }
        )
        .toList()
        .filter { passes -> !passes.isEmpty() }
        .map { passes -> passes.max { pass -> pass[0] }[1] }

    emit:
    combined_protein = IONQUANT.out.combined_protein           // channel: [ val(meta), path(combined_protein.tsv) ]
    combined_peptide = IONQUANT.out.combined_peptide           // channel: [ val(meta), path(combined_peptide.tsv) ]
    combined_ion     = IONQUANT.out.combined_ions              // channel: [ val(meta), path(combined_ion.tsv) ]
    combined_modified_peptide = IONQUANT.out.combined_modified_peptide  // channel: [ val(meta), path(combined_modified_peptide.tsv) ]
    experiment_annotation = IONQUANT.out.experiment_annotation
        .mix(IONQUANT_ISOBARIC.out.experiment_annotation)    // channel: [ val(meta), path(experiment_annotation.tsv) ]
    tmt_abundance    = TMTINTEGRATOR.out.abundance             // channel: [ val(meta), path(abundance_*.tsv) ]
    tmt_config       = TMTINTEGRATOR.out.config                // channel: [ val(meta), path(tmt-integrator-conf.yml) ]
    ionquant_dir     = IONQUANT.out.results_dir                // channel: [ val(meta), path(results_dir) ]
    reprint          = ch_reprint                              // channel: [ val(mode), path(reprint.<mode>.tsv) ] - REPRINT tables IonQuant leaves in the work dir
    ionquant_ms1_dir = IONQUANT_MS1.out.results_dir            // channel: [ val(meta), path(results_dir) ]
    tmt_dir          = TMTINTEGRATOR.out.results_dir           // channel: [ val(meta), path(tmt-report) ]
    tmt_unmod_dir    = TMTINTEGRATOR_UNMOD.out.results_dir     // channel: [ val(meta), path(tmt-report-unmod) ]
    isobaric_dir     = IONQUANT_ISOBARIC.out.results_dir       // channel: [ val(meta), path(results_dir) ]
    quant_dir        = ch_quant_dir                            // channel: path(results_dir) - the last IonQuant pass's work dir
}
