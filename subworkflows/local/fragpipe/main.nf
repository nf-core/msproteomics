/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UNIFIED FRAGPIPE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A dynamic nf-core compliant subworkflow that runs FragPipe tools based on
    workflow configuration.

    This subworkflow composes focused subworkflows:
    1. FRAGPIPE_TOOL_CONFIGS - Tool configs from the workflow file (parser) or FragPipe's plan (--fragpipe_planner)
    2. FRAGPIPE_DATABASE - Database preparation (or use prebuilt)
    3. FRAGPIPE_CONVERT - Optional DIA preprocessing (DIA-Umpire, diaTracer)
    4. FRAGPIPE_SEARCH - MSFragger search with native calibration
    5. FRAGPIPE_VALIDATE - Percolator / PeptideProphet / PTMProphet
    6. FRAGPIPE_INFERENCE - ProteinProphet + Philosopher Filter (+ iProphet, Abacus)
    7. FRAGPIPE_PHILOSOPHER_QUANT - FreeQuant / LabelQuant + final Report, ExportMatchedFragments
       then WRITESUBMZML (tab-run.write_sub_mzml) per LC-MS file
    8. FRAGPIPE_GLYCO - PairScans, O-Pair, PTM-Shepherd, MBG (optional, in FragPipe's order)
    9. FRAGPIPE_QUANT - IonQuant / TMT-Integrator
    10. FRAGPIPE_SPECLIB - Spectral library generation + DIA-NN quantification (optional)
    11. FRAGPIPE_EXPORT - SAINT / FPOP / FPOP coADAPTr / Metaproteomics (optional)

    Tool Execution Order (DAG):
    DATABASE + CONVERT (parallel)
        -> SEARCH (MSFragger + optional CrystalC/MSBooster)
        -> VALIDATE (Percolator or PeptideProphet, optional PTMProphet)
        -> INFERENCE (ProteinProphet + Philosopher Filter)
        -> PHILOSOPHER_QUANT (FreeQuant / LabelQuant + Report, ExportMatchedFragments)
        -> GLYCO (PairScans -> O-Pair -> PTM-Shepherd -> MBG; each rewrites psm.tsv)
        -> QUANT / SPECLIB on the rewritten group directories; EXPORT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Tool configs: workflow parser or FragPipe's exported plan (--fragpipe_planner)
include { FRAGPIPE_TOOL_CONFIGS    } from '../fragpipe_tool_configs/main'
include { CHECKCENTROID           } from '../../../modules/local/checkcentroid/main'

// Import utility functions
include { shouldRunTool; requireToolField; msboosterFraggerParamsType; checkTimsTofTools; fragpipeWorkflowProperties; checkCentroidFile; getMsfraggerNumSlices; checkFpopTools } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

// Import subworkflows
include { FRAGPIPE_DATABASE                                    } from '../fragpipe_database/main'
include { FRAGPIPE_CONVERT                                     } from '../fragpipe_convert/main'
include { FRAGPIPE_SEARCH                                      } from '../fragpipe_search/main'
include { FRAGPIPE_VALIDATE  } from '../fragpipe_validate/main'
include { FRAGPIPE_INFERENCE } from '../fragpipe_inference/main'
include { FRAGPIPE_PHILOSOPHER_QUANT } from '../fragpipe_philosopher_quant/main'
include { FRAGPIPE_QUANT     } from '../fragpipe_quant/main'
include { FRAGPIPE_SPECLIB   } from '../fragpipe_speclib/main'
include { FRAGPIPE_GLYCO     } from '../fragpipe_glyco/main'
include { FRAGPIPE_EXPORT    } from '../fragpipe_export/main'
include { FRAGPIPE_MATCHED_FRAGMENTS_INPUTS } from '../fragpipe_matched_fragments/main'
include { WRITESUBMZML     } from '../../../modules/local/writesubmzml/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FRAGPIPE {
    take:
    ch_file            // channel: [ val(sample_meta), path(mzml/raw/d) ] - sample files (meta.data_type: FragPipe data type)
    ch_database        // channel: [ val(meta), path(fasta) ] - FASTA to update with decoys (empty unless --add_philosopher_decoy / --add_philosopher_contaminants)
    ch_workflow_file   // channel: [ val(meta), path(workflow_file) ] - FragPipe workflow file
    ch_decoy_db        // channel: [ val(meta), path(fasta) ] - FASTA searched as given (empty when ch_database is used)
    ch_annotation      // channel: path(annotation_file) - value channel, [] when not used
    ch_manifest        // channel: path(fragpipe-files.fp-manifest) - value channel (fragpipeManifestEntries)
    ch_taxon_files     // channel: value [ path(names.dmp) or [], path(nodes.dmp) or [] ] - Metaproteomics taxonomy dumps; [] = FragPipe's tools/metaproteomics/{names,nodes}.dmp (CmdMetaproteomics.java:39-42)

    main:

    //
    // STEP 1: Parse workflow config to get tool run flags and args
    //
    FRAGPIPE_TOOL_CONFIGS(ch_file, ch_workflow_file, ch_decoy_db, ch_annotation, params.fragpipe_planner)

    // Read tool configs JSON as a string value channel.
    // Strings avoid the Groovy Map-as-named-parameter bug in channel operators.
    // Use splitText to avoid blocking I/O inside map{} on cloud executors.
    // FragPipe aborts when Crystal-C, FreeQuant or LabelQuant would run on timsTOF data
    // (FragpipeRun.java:2420-2443, see checkTimsTofTools).
    ch_tool_configs = FRAGPIPE_TOOL_CONFIGS.out.tool_configs_json
        .map { _meta, json_file -> json_file }
        .splitText()
        .collect()
        .map { lines -> lines.join('') }
        .combine(ch_file.map { _meta, spectra -> spectra.name }.toList().map { names -> [names] })
        .combine(ch_workflow_file.map { _meta, workflow_file -> fragpipeWorkflowProperties(workflow_file).getProperty('workflow.input.data-type.im-ms') == 'true' })
        .map { configs, names, im_ms -> checkTimsTofTools(configs, names, im_ms) }
        .first()

    // Extract config files from FRAGPIPE_TOOL_CONFIGS as value channels.
    // These are params files for tools that use them (MSFragger, MSBooster, etc.).
    // The parser writes one config per tool prefix present in the workflow file, and the four
    // MSFragger params files always; the planner writes the params files FragPipe wrote.
    // MSFragger params per data type: [[params_type, path], ...] (CmdMsfragger.java:444-457); the parser
    // writes all four, with FragPipe's closed-search defaults when the workflow has no msfragger.* keys;
    // the planner one per data type of the manifest.
    ch_fragger_params = FRAGPIPE_TOOL_CONFIGS.out.msfragger_config.map { _meta, params -> ['DDA', params] }
        .mix(FRAGPIPE_TOOL_CONFIGS.out.msfragger_dia_config.map { _meta, params -> ['DIA', params] })
        .mix(FRAGPIPE_TOOL_CONFIGS.out.msfragger_gpfdia_config.map { _meta, params -> ['GPF-DIA', params] })
        .mix(FRAGPIPE_TOOL_CONFIGS.out.msfragger_dda_plus_config.map { _meta, params -> ['DDA+', params] })
        .toList()
        // empty only when PARSEFRAGPIPEWORKFLOW failed: the run then stops with the parser's error, which
        // a NullPointerException on sets[0] downstream would hide
        .filter { sets -> !sets.isEmpty() }

    ch_msbooster_params = FRAGPIPE_TOOL_CONFIGS.out.msbooster_config
        .map { _meta, params -> params }
        .first()


    // MSBooster reads one MSFragger params file chosen from the data types of ALL manifest files
    // (TabWorkflow.hasDataType) and whether DIA-Umpire / diaTracer run (CmdMSBooster.java:134-149).
    ch_msbooster_fragger_params = ch_file
        .map { meta, _spectra -> meta.data_type }
        .toList()
        .map { data_types -> [data_types] }
        .combine(ch_tool_configs)
        .combine(ch_fragger_params.map { sets -> [sets] })
        .map { data_types, configs, sets ->
            if (!shouldRunTool(configs, 'msbooster')) {
                return []
            }
            def params_type = msboosterFraggerParamsType(data_types, shouldRunTool(configs, 'diaumpire'), shouldRunTool(configs, 'diatracer'))
            def fragger_params = sets.find { set -> set[0] == params_type }
            if (!fragger_params) {
                error("MSBooster reads the MSFragger params of data type ${params_type}, but the parser wrote only ${sets.collect { set -> set[0] }}")
            }
            fragger_params[1]
        }
        .first()  // value channel: combine() above emits a queue channel

    // FragPipe leaves DIA-Quant files out of every stage from DIA-Umpire to SpecLibGen
    // (FragpipeRun.java:1273-1280); only DIA-NN (FRAGPIPE_SPECLIB) quantifies them.
    //
    // FragPipe's CheckCentroid on the first LC-MS file gates DIA preprocessing and the search
    // (FragpipeRun.java:1333-1340; graph :2448-2451)
    //
    // Skipped when no file is searched (only DIA-Quant files), as FragPipe (FragpipeRun.java:1335)
    CHECKCENTROID(
        ch_file.toList()
            .flatMap { files ->
                def centroid_input = checkCentroidFile(files)
                centroid_input != null ? [centroid_input] : []
            }
    )

    ch_fragpipe_files = ch_file
        .filter { meta, _spectra -> meta.data_type != 'DIA-Quant' }
        .combine(CHECKCENTROID.out.log.map { _meta, log -> log })
        .map { meta, spectra, _log -> [meta, spectra] }

    //
    // STEP 2: FragPipe's database decoy check (CHECKDBCONFIG); ch_database only holds a
    // FASTA when --add_philosopher_decoy / --add_philosopher_contaminants asks for FragPipe's database update first.
    //
    FRAGPIPE_DATABASE(ch_database, ch_decoy_db, ch_tool_configs)
    ch_fasta = FRAGPIPE_DATABASE.out.fasta_path

    //
    // STEP 3: Optional DIA preprocessing (per-sample, parallel); LC-MS files are searched as given
    //
    FRAGPIPE_CONVERT(ch_fragpipe_files, ch_tool_configs)

    ch_mzml = FRAGPIPE_CONVERT.out.mzml

    //
    // STEP 4: Database search (MSFragger with native calibration)
    //
    FRAGPIPE_SEARCH(
        ch_mzml,
        ch_fasta,
        ch_tool_configs,
        ch_fragger_params,
        ch_msbooster_fragger_params,
        ch_msbooster_params,
        ch_tool_configs.map { configs -> getMsfraggerNumSlices(configs) }
    )
    ch_search_pepxml       = FRAGPIPE_SEARCH.out.pepxml
    ch_search_pin          = FRAGPIPE_SEARCH.out.pin
    ch_search_pin_rescored = FRAGPIPE_SEARCH.out.pin_rescored

    // FragPipe's LC-MS directory after the search (LC-MS file + MSFragger's mzML/mzBIN outputs),
    // one [meta, file] item per file, for the tools that read spectra directories;
    // Percolator/PeptideProphet/glyco name one file per run.
    ch_mzml_downstream = FRAGPIPE_SEARCH.out.lcms
        .flatMap { meta, files -> files.collect { f -> [meta, f] } }
    ch_searched_lcms = FRAGPIPE_SEARCH.out.lcms.map { meta, files -> [meta, files[0]] }

    //
    // STEP 5: Percolator / PeptideProphet / PTMProphet validation
    //
    FRAGPIPE_VALIDATE(
        ch_search_pin,                        // MSFragger PIN (PercolatorOutputToPepXML reads it)
        ch_search_pin_rescored,               // MSBooster _edited.pin, or the MSFragger PIN
        ch_search_pepxml,                     // pepXML (Crystal-C-filtered or raw)
        ch_searched_lcms,                     // searched LC-MS file
        FRAGPIPE_SEARCH.out.calibrated_mzml,
        FRAGPIPE_SEARCH.out.uncalibrated_mzml,
        ch_fasta,
        ch_tool_configs
    )

    //
    // STEP 6: ProteinProphet + Philosopher Filter (aggregate -> per-sample)
    //
    def aggregate_meta = [id: 'fragpipe']

    FRAGPIPE_INFERENCE(
        FRAGPIPE_VALIDATE.out.pepxml,
        FRAGPIPE_VALIDATE.out.pepxml_group,
        ch_fasta,
        ch_tool_configs,
        aggregate_meta
    )

    //
    // STEP 6b: FreeQuant / LabelQuant + final Philosopher Report, then Export matched fragments
    // (FragpipeRun.java:2463-2471, :2490). Every later tool reads these tables: the glyco chain
    // (IN-edges from PhilosopherReport), IonQuant, TMT-Integrator, SpecLibGen and the exports.
    //
    // Spectra the export reads: the LC-MS directory after the search (searched mzML files and
    // MSFragger's _uncalibrated.mzML of .raw / .d files)
    FRAGPIPE_MATCHED_FRAGMENTS_INPUTS(
        ch_file,
        ch_mzml_downstream,
        FRAGPIPE_CONVERT.out.umpire_mzml
    )

    FRAGPIPE_PHILOSOPHER_QUANT(
        FRAGPIPE_INFERENCE.out.results_dir,
        ch_mzml_downstream,
        ch_file,
        ch_annotation,
        ch_tool_configs,
        FRAGPIPE_MATCHED_FRAGMENTS_INPUTS.out.inputs,
        aggregate_meta
    )
    ch_results_dir = FRAGPIPE_PHILOSOPHER_QUANT.out.final_results_dir

    //
    // STEP 6c: Write sub mzML (tab-run.write_sub_mzml), per LC-MS file after the final Philosopher
    // report, like ExportMatchedFragments (FragpipeRun.java:2388-2404, graph :2489). WriteSubMzml
    // reads the file as given, MSFragger's <stem>.mzBIN_all next to it and its group's psm.tsv
    // (CmdWriteSubMzml.java:97-130); it uses the Spectrum and Probability columns only, which the
    // export and the later psm.tsv rewrites keep.
    //
    ch_group_psm = ch_results_dir
        .map { meta, dir -> [meta.group, dir.resolve('psm.tsv')] }

    ch_for_write_sub_mzml = ch_fragpipe_files
        .map { meta, spectra -> [meta.id, meta, spectra] }
        .join(FRAGPIPE_SEARCH.out.mzbin_all.map { meta, mzbin -> [meta.id, mzbin] }, failOnDuplicate: true, remainder: true)
        .combine(ch_tool_configs)
        .filter { _key, _meta, _spectra, _mzbin, configs -> shouldRunTool(configs, 'writesubmzml') }
        .map { key, meta, spectra, mzbin, configs ->
            if (mzbin == null) {
                error("WRITESUBMZML: MSFragger wrote no .mzBIN_all for ${key}, but Write sub mzML is on (write_mzbin_all = 1, FragpipeRun.java:1481)")
            }
            [meta.group, meta, spectra, mzbin, requireToolField(configs, 'writesubmzml', 'prob_threshold')]
        }
        // every searched file's group has a Report psm.tsv; combine(by:) would drop a file without one silently
        .combine(ch_group_psm.map { group, _psm -> group }.collect().ifEmpty([]).map { groups -> [groups] })
        .map { group, meta, spectra, mzbin, threshold, groups ->
            if (!(group in groups)) {
                error("WRITESUBMZML: group '${group}' of ${meta.id} has no Report psm.tsv (groups with one: ${groups})")
            }
            [group, meta, spectra, mzbin, threshold]
        }
        .combine(ch_group_psm, by: 0)
        .map { _group, meta, spectra, mzbin, threshold, psm -> [meta, spectra, mzbin, psm, threshold] }

    WRITESUBMZML(ch_for_write_sub_mzml)

    //
    // STEP 7a: Glyco chain - PairScans -> O-Pair -> PTM-Shepherd -> MBG (FragpipeRun.java:2472-2477),
    // on the final per-group tables (after the final Report and ExportMatchedFragments). O-Pair, PTM-Shepherd and MBG
    // rewrite psm.tsv; IonQuant and SpecLibGen run after them on the rewritten group directories
    // (FragpipeRun.java:2477, 2483).
    //
    FRAGPIPE_GLYCO(
        ch_results_dir,
        ch_fragpipe_files,
        ch_searched_lcms,                     // one searched LC-MS file per run (PairScans, O-Pair, PTM-Shepherd, MBG)
        FRAGPIPE_SEARCH.out.calibrated_mzml,
        FRAGPIPE_SEARCH.out.uncalibrated_mzml,
        ch_fasta,
        ch_tool_configs,
        aggregate_meta
    )
    ch_group_dirs = FRAGPIPE_GLYCO.out.psm_dirs

    // AppendFile input for IonQuant's modlist ([] when AppendFile does not run)
    ch_glyco_masses = FRAGPIPE_GLYCO.out.modlist_append

    //
    // STEP 7b: Quantification - IonQuant / TMT-Integrator (aggregate) on the tables the glyco chain
    // left; every later reader uses them, as FragPipe's tools read the work dir's psm.tsv.
    //
    FRAGPIPE_QUANT(
        ch_group_dirs,
        ch_mzml_downstream,
        ch_file,
        ch_annotation,
        ch_tool_configs,
        ch_glyco_masses,
        aggregate_meta
    )

    ch_final_results = ch_group_dirs
    ch_psms = ch_final_results
        .map { meta, dir -> [meta, dir.resolve('psm.tsv')] }
    ch_peptides = ch_final_results
        .map { meta, dir -> [meta, dir.resolve('peptide.tsv')] }
    ch_proteins = ch_final_results
        .map { meta, dir -> [meta, dir.resolve('protein.tsv')] }

    //
    // STEP 7c: Spectral library generation and DIA-NN quantification (optional, aggregate), on the
    // psm.tsv the glyco chain left
    //
    FRAGPIPE_SPECLIB(
        ch_psms,
        ch_peptides,
        ch_proteins,
        FRAGPIPE_VALIDATE.out.pepxml,
        ch_searched_lcms,                     // files MSFragger searched
        FRAGPIPE_SEARCH.out.uncalibrated_mzml,
        ch_file,
        ch_manifest,
        ch_fasta,
        ch_tool_configs,
        aggregate_meta
    )

    //
    // STEP 7e: Export - SAINT, FPOP, FPOP coADAPTr, Metaproteomics (optional)
    //
    // FragPipe's FPOP checks that need the LC-MS data types (checkFpopTools)
    ch_export_configs = ch_file
        .map { meta, _spectra -> meta.data_type }
        .toList()
        .map { data_types -> [data_types] }
        .combine(ch_tool_configs)
        .map { data_types, configs -> checkFpopTools(configs, data_types) }
        .first()

    // Tables FPOP coADAPTr edits, by fpopCoadaptrSource: [source, layout, [dirs]]. FragPipe walks
    // its work dir; here the IonQuant results directory is passed whole (every psm.tsv in it),
    // the Philosopher results per group (their psm.tsv).
    ch_psm_tables = FRAGPIPE_QUANT.out.ionquant_dir
        .map { _meta, dir -> ['ionquant', 'tree', [dir]] }
        .mix(
            FRAGPIPE_QUANT.out.isobaric_dir.map { _meta, dir -> ['isobaric', 'tree', [dir]] },
            ch_final_results
                .filter { meta, _dir -> meta.report_after_quant }
                .map { _meta, dir -> dir }
                .collect(sort: true)
                .map { dirs -> ['report', 'groups', dirs] },
            ch_final_results
                .filter { meta, _dir -> !meta.report_after_quant }
                .map { _meta, dir -> dir }
                .collect(sort: true)
                .map { dirs -> ['filter', 'groups', dirs] }
        )

    // DIA-NN report for FPOP coADAPTr's DIA mode (FPOP_DiannReport_Editor edits
    // dia-quant-output/report.tsv): [ val(meta), path(report.tsv) ], empty when DIA-NN does not run
    ch_diann_report = FRAGPIPE_SPECLIB.out.diann_report

    // FPOP Quant input (CmdFpopQuant.java:56-63): IonQuant's combined_modified_peptide.tsv,
    // or in TMT mode the tmt-report and tmt-report-unmod TMT-Integrator directories.
    ch_fpop_inputs = FRAGPIPE_QUANT.out.combined_modified_peptide
        .map { _meta, modpep -> [modpep, [], []] }
        .mix(
            FRAGPIPE_QUANT.out.tmt_dir
                .combine(ch_export_configs)
                .filter { _meta, _tmt_report, configs -> shouldRunTool(configs, 'tmtintegrator_unmod') }
                .map { meta, tmt_report, _configs -> [meta, tmt_report] }
                .join(FRAGPIPE_QUANT.out.tmt_unmod_dir, failOnMismatch: true, failOnDuplicate: true)
                .map { _meta, tmt_report, tmt_report_unmod -> [[], tmt_report, tmt_report_unmod] }
        )

    // Metaproteomics project dir: FP-Meta reads every <run>.pin (MSFragger) and
    // <run>_percolator_{target,decoy}_psms.tsv of the work dir, right after Percolator
    // (FragpipeRun.java:2458).
    ch_metaproteomics_files = ch_search_pin
        .map { _meta, pin -> pin }
        .mix(
            FRAGPIPE_VALIDATE.out.target_psms.map { _meta, psms -> psms },
            FRAGPIPE_VALIDATE.out.decoy_psms.map { _meta, psms -> psms }
        )
        .collect(sort: true)

    // SAINTexpress reads reprint.<mode>.tsv from the work dir (CmdSaintExpress.java:59): written by
    // Abacus --reprint, or by IonQuant when it runs (Abacus then does not, FragpipeRun.java:1786-1788)
    ch_saint_files = FRAGPIPE_INFERENCE.out.reprint_spc.map { _meta, tsv -> ['spc', tsv] }
        .mix(FRAGPIPE_INFERENCE.out.reprint_int.map { _meta, tsv -> ['int', tsv] })
        .mix(FRAGPIPE_QUANT.out.reprint)

    FRAGPIPE_EXPORT(
        ch_fasta,
        ch_psm_tables,
        ch_diann_report,
        ch_fpop_inputs,
        ch_metaproteomics_files,
        ch_taxon_files,
        ch_saint_files,
        ch_export_configs,
        aggregate_meta
    )

    emit:
    // LC-MS files with MSFragger's side files
    searched_lcms         = ch_searched_lcms                                   // channel: [ val(meta), path(file) ] - the LC-MS file MSFragger searched
    mzml                  = ch_mzml_downstream                                 // channel: [ val(meta), path(file) ] - LC-MS files and MSFragger's side files

    // Search outputs
    pepxml                = ch_search_pepxml                                   // channel: [ val(meta), path(pepxml) ]
    pin                   = ch_search_pin                                      // channel: [ val(meta), path(pin) ]
    database              = ch_fasta                                           // channel: path(fasta) - the run database (decoys, contaminants)
    msfragger_logs        = FRAGPIPE_SEARCH.out.msfragger_logs                 // channel: [ val(meta), path(log) ] - MSFragger console output
    sub_mzml              = WRITESUBMZML.out.sub_mzml                        // channel: [ val(meta), path(<stem>_sub.mzML) ] - Write sub mzML

    // Validation outputs
    validated_pepxml      = FRAGPIPE_VALIDATE.out.pepxml                       // channel: [ val(meta), path(pepxml) ]

    // Inference outputs
    protxml               = FRAGPIPE_INFERENCE.out.protxml                     // channel: [ val(meta), path(combined.prot.xml) ]
    psms                  = ch_psms                                            // channel: [ val(meta), path(psm.tsv) ] - with matched fragments when exported
    peptides              = ch_peptides                                        // channel: [ val(meta), path(peptide.tsv) ]
    proteins              = ch_proteins                                        // channel: [ val(meta), path(protein.tsv) ]
    results_dir           = ch_final_results                                   // channel: [ val(meta), path(results_dir) ]

    // Quantification outputs
    combined_protein      = FRAGPIPE_QUANT.out.combined_protein                // channel: [ val(meta), path(combined_protein.tsv) ]
    combined_peptide      = FRAGPIPE_QUANT.out.combined_peptide                // channel: [ val(meta), path(combined_peptide.tsv) ]
    combined_ion          = FRAGPIPE_QUANT.out.combined_ion                    // channel: [ val(meta), path(combined_ion.tsv) ]
    ionquant_dir          = FRAGPIPE_QUANT.out.ionquant_dir                    // channel: [ val(meta), path(results_dir) ]
    tmt_abundance         = FRAGPIPE_QUANT.out.tmt_abundance                   // channel: [ val(meta), path(abundance_*.tsv) ]
    tmt_dir               = FRAGPIPE_QUANT.out.tmt_dir                         // channel: [ val(meta), path(tmt-report) ]
    combined_modified_peptide = FRAGPIPE_QUANT.out.combined_modified_peptide   // channel: [ val(meta), path(combined_modified_peptide.tsv) ]

    // PTM / glyco outputs
    ptm_results           = FRAGPIPE_GLYCO.out.ptm_results                     // channel: [ val(meta), path(ptm-shepherd-output) ]
    glyco_psm_dirs        = FRAGPIPE_GLYCO.out.psm_dirs                        // channel: [ val(meta), path(group_dir) ]

    // Spectral library and DIA-NN outputs
    speclib               = FRAGPIPE_SPECLIB.out.library_tsv                   // channel: [ val(meta), path(library.tsv) ]
    diann_report          = FRAGPIPE_SPECLIB.out.diann_report                  // channel: [ val(meta), path(report.tsv) ]
    diann_pg_matrix       = FRAGPIPE_SPECLIB.out.diann_pg_matrix               // channel: [ val(meta), path(report.pg_matrix.tsv) ]
    diann_results         = FRAGPIPE_SPECLIB.out.diann_results                 // channel: [ val(meta), [ path(dia-quant-output/*) ] ]

    // Export outputs
    taxonomy_results      = FRAGPIPE_EXPORT.out.taxonomy_results               // channel: [ val(meta), path(results_dir) ]
    fpop_results          = FRAGPIPE_EXPORT.out.fpop_results                   // channel: [ val(meta), path(results_dir) ]
    fpop_psms             = FRAGPIPE_EXPORT.out.fpop_psms                      // channel: [ val(meta), path(psm.tsv) ]
    fpop_diann_report     = FRAGPIPE_EXPORT.out.fpop_diann_report              // channel: [ val(meta), path(report.tsv) ]

}
