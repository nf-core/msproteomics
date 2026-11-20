/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_GLYCO SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The glyco chain between Philosopher Report and IonQuant. FragPipe's graph edges
    (FragpipeRun.java:2472-2477, 24.0):

        addToGraph(graphOrder, cmdPairScans, DIRECTION.IN, cmdPhilosopherReport, cmdPhilosopherAbacus);
        addToGraph(graphOrder, cmdOPair, DIRECTION.IN, cmdPairScans);
        addToGraph(graphOrder, cmdPtmshepherd, DIRECTION.IN, cmdPhilosopherReport, cmdPhilosopherAbacus);
        addToGraph(graphOrder, cmdAppendFile, DIRECTION.IN, cmdPtmshepherd);
        addToGraph(graphOrder, cmdMBGMatch, DIRECTION.IN, cmdPhilosopherReport, cmdPhilosopherAbacus, cmdPtmshepherd, cmdOPair);
        addToGraph(graphOrder, cmdIonquant, DIRECTION.IN, cmdPhilosopherReport, cmdPhilosopherAbacus, cmdPtmshepherd, cmdMBGMatch);

    PTM-Shepherd does not depend on O-Pair, but FragPipe never runs both (FragpipeRun.java:2099-2111,
    `if (oPairPanel.isRun() && isRunShepherd) ... return false;`, so FragPipe's plan stops such a workflow).
    So at most one of them rewrites the Report psm.tsv, and the chain below gives every tool the
    psm.tsv it reads in FragPipe:

        PairScans -> O-Pair  |  PTM-Shepherd (+ AppendFile masses)   ->  MBG  ->  IonQuant / SpecLibGen

    O-Pair, PTM-Shepherd and MBG each rewrite the group's psm.tsv in place in FragPipe, so here every
    stage works on a copy of the group directory and the next stage reads that copy. A stage that
    does not run passes the group directories on unchanged. The emitted psm_dirs are what IonQuant
    reads (declared edge, L2477). SpecLibGen's declared edges are only Report and O-Pair (L2483:
    `addToGraph(graphOrder, cmdSpecLibGen, DIRECTION.IN, cmdPhilosopherReport, cmdOPair)`); that it
    also runs after PTM-Shepherd and MBG is the order FragPipe's topological iteration gives, inferred
    and still to be confirmed with a FragPipe --dry-run "Execution order" of a glyco-*-DIA workflow.

    Modules:
    - PAIRSCANS (per LC-MS file), OPAIR (per group), PTMSHEPHERD (all groups), MBG (per group)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PAIRSCANS    } from '../../../modules/local/pairscans/main'
include { OPAIR        } from '../../../modules/local/opair/main'
include { MBG          } from '../../../modules/local/mbg/main'
include { FRAGPIPE_PTM } from '../fragpipe_ptm/main'
include { shouldRunTool; getToolArgs; getToolField; getToolEntry; mbgNoPasef; opairSpectraSource } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

workflow FRAGPIPE_GLYCO {
    take:
    ch_results_dir        // channel: [ val(meta), path(group_dir) ] - per FragPipe group (meta.group), psm.tsv inside
    ch_lcms               // channel: [ val(meta), path(spectra) ] - LC-MS files in FragPipe file-table order (meta.id = basename)
    ch_mzml               // channel: [ val(meta), path(spectra) ] - the LC-MS files MSFragger searched (meta.group)
    ch_calibrated_mzml    // channel: [ val(meta), path(*_calibrated.mzML) ] - MSFragger write_calibrated_mzml
    ch_uncalibrated_mzml  // channel: [ val(meta), path(*_uncalibrated.mzML) ] - MSFragger .raw/.d
    ch_fasta              // channel: path(fasta) - the searched database
    ch_tool_configs       // channel: val(tool_configs_json)
    aggregate_meta        // val(meta) for aggregate outputs

    main:

    //
    // PairScans: <base>.pairs per LC-MS file (CmdPairScans.java:108-128)
    //
    PAIRSCANS(
        ch_mzml
            .combine(ch_tool_configs)
            .filter { _meta, _spectra, configs -> shouldRunTool(configs, 'pairscans') }
            .map { meta, spectra, configs -> [meta, spectra, getToolArgs(configs, 'pairscans')] }
    )

    //
    // O-Pair: per group, on the spectra its filelist names (CmdOPair.java:99-117)
    //
    ch_opair_spectra = ch_mzml
        .combine(ch_tool_configs)
        .filter { _meta, _spectra, configs -> shouldRunTool(configs, 'opair') }
        .map { meta, spectra, configs ->
            [meta.id, meta, spectra, opairSpectraSource(getToolField(configs, 'opair', 'spectra'), spectra.name)]
        }
        .branch { _id, _meta, _spectra, source ->
            calibrated: source == 'calibrated'
            uncalibrated: source == 'uncalibrated'
            input: true
        }

    // remainder: MSFragger also writes (un)calibrated mzML that O-Pair does not list; those rows
    // have no left side and are dropped. A listed file MSFragger did not write is an error.
    ch_opair_written = ch_opair_spectra.calibrated
        .join(ch_calibrated_mzml.map { meta, mzml -> [meta.id, mzml] }, failOnDuplicate: true, remainder: true)
        .mix(
            ch_opair_spectra.uncalibrated
                .join(ch_uncalibrated_mzml.map { meta, mzml -> [meta.id, mzml] }, failOnDuplicate: true, remainder: true)
        )
        .filter { item -> item[1] != null }
        .map { _id, meta, spectra, source, listed ->
            if (listed == null) {
                error("O-Pair lists the ${source} mzML of ${spectra.name}, but MSFragger did not write one")
            }
            [meta, listed]
        }

    // FragPipe lists every LC-MS file of the group, sorted (CmdOPair.java:100)
    ch_opair_group_files = ch_opair_written
        .mix(ch_opair_spectra.input.map { _id, meta, spectra, _source -> [meta, spectra] })
        .map { meta, listed -> [meta.group, listed] }
        .groupTuple()
        .map { group, listed -> [group, listed.sort { f -> f.name }] }

    ch_group_pairs = PAIRSCANS.out.pairs
        .map { meta, pairs -> [meta.group, pairs] }
        // by file name: a fixed order for the OPAIR task hash (resume); the files sit in different task dirs
        .groupTuple(sort: { a, b -> a.name <=> b.name })

    // the user's oxonium ion definitions file (opair.oxonium_filtering_file), [] for FragPipe's default in the tools dir;
    // a declared file that does not exist stops the run; wrapped in a list so combine adds [] as one element
    ch_oxonium_file = ch_tool_configs.map { configs ->
        def oxonium_file = getToolEntry(configs, 'opair').oxonium_file
        [oxonium_file ? file(oxonium_file, checkIfExists: true) : []]
    }

    ch_for_opair = ch_results_dir
        .combine(ch_tool_configs)
        .combine(ch_oxonium_file)
        .filter { _meta, _dir, configs, _oxonium_file -> shouldRunTool(configs, 'opair') }
        .map { meta, dir, configs, oxonium_file ->
            [meta.group, meta, dir, getToolArgs(configs, 'opair'), oxonium_file]
        }
        .join(ch_opair_group_files, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_group_pairs, remainder: true)
        .filter { item -> item[1] != null }
        .map { _group, meta, dir, argv, oxonium_file, listed, pairs ->
            [meta, dir, listed, pairs ?: [], oxonium_file, listed.collect { f -> f.name }, argv]
        }

    OPAIR(ch_for_opair)

    ch_after_opair = OPAIR.out.psm_dir
        .mix(
            ch_results_dir
                .combine(ch_tool_configs)
                .filter { _meta, _dir, configs -> !shouldRunTool(configs, 'opair') }
                .map { meta, dir, _configs -> [meta, dir] }
        )

    //
    // PTM-Shepherd (profiling and/or glycan assignment) after O-Pair
    //
    FRAGPIPE_PTM(
        ch_after_opair,
        ch_mzml,
        ch_tool_configs,
        aggregate_meta,
        ch_fasta
    )

    //
    // MBG: per group, after PTM-Shepherd and O-Pair (CmdMBGMatch.java:125-143); reads every
    // LC-MS file of the run through fragpipe-files.fp-manifest
    //
    ch_lcms_order = ch_lcms
        .map { meta, spectra -> [meta.id, spectra.name] }
        .toList()
        .map { order -> [order] }

    ch_mbg_run = ch_mzml
        .toList()
        .map { files -> [files] }
        .combine(ch_lcms_order)
        .combine(ch_tool_configs)
        .filter { _files, _order, configs -> shouldRunTool(configs, 'mbg') }
        .map { files, order, configs ->
            def by_id = files.collectEntries { item -> [item[0].id, item] }
            def ordered = order.collect { id_name ->
                if (!by_id.containsKey(id_name[0])) {
                    error("MBG manifest: LC-MS file ${id_name[1]} has no searched spectra")
                }
                by_id[id_name[0]]
            }
            def rows = ordered.collect { item -> [item[1].name, item[0].experiment, item[0].bioreplicate, item[0].data_type] }
            def nopasef = mbgNoPasef(order.collect { id_name -> id_name[1] }, getToolField(configs, 'mbg', 'im_ms', false))
            def argv = getToolArgs(configs, 'mbg').collect { token -> token == '@NOPASEF@' ? nopasef : token }
            [ordered.collect { item -> item[1] }, rows, argv]
        }

    MBG(
        FRAGPIPE_PTM.out.psm_dirs
            .combine(ch_mbg_run)
            .map { meta, dir, spectra, rows, argv -> [meta, dir, spectra, rows, argv] }
    )

    ch_psm_dirs = MBG.out.psm_dir
        .mix(
            FRAGPIPE_PTM.out.psm_dirs
                .combine(ch_tool_configs)
                .filter { _meta, _dir, configs -> !shouldRunTool(configs, 'mbg') }
                .map { meta, dir, _configs -> [meta, dir] }
        )

    // AppendFile: PTM-S glyco_masses_list.txt appended to IonQuant's modlist (FragpipeRun.java:1822-1830)
    ch_modlist_append = FRAGPIPE_PTM.out.glyco_masses
        .map { _meta, masses -> masses }
        .toList()
        .map { files -> [files] }
        .combine(ch_tool_configs)
        .map { files, configs ->
            if (!shouldRunTool(configs, 'appendfile')) {
                return []
            }
            if (files.size() != 1) {
                error("AppendFile needs PTM-Shepherd's ptm-shepherd-output/glyco_masses_list.txt, but PTM-Shepherd produced ${files.size()} (glycan assignment must run)")
            }
            files[0]
        }

    emit:
    modlist_append   = ch_modlist_append                 // channel: path(glyco_masses_list.txt) or [] - value channel for IONQUANT
    psm_dirs         = ch_psm_dirs                       // channel: [ val(meta), path(group_dir) ] - per group, after the glyco chain
    ptm_results      = FRAGPIPE_PTM.out.results_dir      // channel: [ val(meta), path(ptm-shepherd-output) ]
    glyco_masses     = FRAGPIPE_PTM.out.glyco_masses     // channel: [ val(meta), path(glyco_masses_list.txt) ]
    shepherd_config  = FRAGPIPE_PTM.out.config           // channel: [ val(meta), path(shepherd.config) ]
    pairs            = PAIRSCANS.out.pairs               // channel: [ val(meta), path(*.pairs) ]
    opair_dirs       = OPAIR.out.psm_dir                 // channel: [ val(meta), path(group_dir) ]
    mbg_dirs         = MBG.out.psm_dir                   // channel: [ val(meta), path(group_dir) ]
    mbg_manifest     = MBG.out.manifest                  // channel: [ val(meta), path(fragpipe-files.fp-manifest) ]
    mbg_log          = MBG.out.log                       // channel: [ val(meta), path(mbg.log) ]
}
