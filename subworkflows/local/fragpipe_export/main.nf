/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_EXPORT SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Export and specialized analysis modules.

    Modules:
    - SAINTEXPRESS (AP-MS scoring for protein-protein interactions)
    - FPOP (Fast Photochemical Oxidation of Proteins analysis)
    - FPOPCOADAPTR (FPOP Modifications column in psm.tsv)
    - METAPROTEOMICS (taxonomic assignment and analysis)

    Execution: Aggregate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REPRINTTOSAINTINPUT } from '../../../modules/local/reprinttosaintinput/main'
include { SAINTEXPRESS   } from '../../../modules/local/saintexpress/main'
include { FPOP           } from '../../../modules/local/fpop/main'
include { FPOPCOADAPTR } from '../../../modules/local/fpopcoadaptr/main'
include { METAPROTEOMICS } from '../../../modules/local/metaproteomics/main'
include { shouldRunTool; getToolArgs; getDecoyTag; requireToolField; fpopCoadaptrSource } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

workflow FRAGPIPE_EXPORT {
    take:
    ch_fasta         // channel: path(fasta) - value channel for database
    ch_psm_tables    // channel: [ val(source), val(layout), [ path(dir), ... ] ] - psm.tsv tables per fpopCoadaptrSource (FPOP coADAPTr)
    ch_diann_report  // channel: [ val(meta), path(dia-quant-output/report.tsv) ] - DIA-NN report (FPOP coADAPTr DIA mode)
    ch_fpop_inputs   // channel: [ path(modpep_tsv) or [], path(tmt-report) or [], path(tmt-report-unmod) or [] ] - FPOP Quant
    ch_meta_files    // channel: [ path(file), ... ] - MSFragger .pin + Percolator *_psms.tsv (Metaproteomics)
    ch_taxon_files   // channel: [ path(names.dmp) or [], path(nodes.dmp) or [] ] - Metaproteomics taxonomy
    ch_saint_files   // channel: [ val(mode), path(reprint.<mode>.tsv) ] - REPRINT tables of the work dir, for SAINT
    ch_tool_configs  // channel: val(tool_configs_map) - parsed JSON config
    aggregate_meta   // val(meta) for aggregate outputs

    main:

    // Decoy tag (workflow database.decoy-tag) - single source for every decoy-aware tool
    ch_decoy_tag = ch_tool_configs.map { configs -> getDecoyTag(configs) }

    //
    // SAINTEXPRESS: FragPipe's Downstream run over both modes (TabDownstream.java:137-145) when the workflow
    // sets saintexpress.run-saint-express=true (a deliberate difference, see FRAGPIPE_INFERENCE STEP 3). Each mode reads the
    // work dir's reprint.<mode>.tsv, whoever wrote it (Abacus --reprint or IonQuant). A missing table makes
    // CmdSaintExpress.configure log an error and return false, and only that mode is skipped
    // (CmdSaintExpress.java:72-79):
    //     log.error("Input file " + inputPath.toAbsolutePath() + " is not accessible.");
    // REPRINTTOSAINTINPUT writes the inter/bait/prey.dat inputs (CmdSaintExpress.reprintTsvToSaintInput), or the
    // message (task dir = work dir) to saintexpress-<mode>/skipped_mode.txt; it is also logged here.
    //
    ch_saint_on = ch_tool_configs.filter { configs -> shouldRunTool(configs, 'saintexpress') }
    ch_saint_on.subscribe { _configs ->
        log.info("SAINTexpress runs as an optional post-run step: FragPipe's Downstream-tab action, not part of FragPipe's run (results in downstream/saintexpress/)")
    }

    ch_for_saintexpress = ch_saint_files
        .toList()
        .map { items -> [items] }
        .combine(ch_saint_on)
        .flatMap { items, configs ->
            ['spc', 'int'].collect { mode ->
                def found = items.findAll { item -> item[0] == mode }
                if (found.size() > 1) {
                    error("SAINTexpress: more than one reprint.${mode}.tsv for the work dir: ${found.collect { item -> item[1] }}")
                }
                if (found && found[0][1].name != "reprint.${mode}.tsv") {
                    error("SAINTexpress: ${mode} mode reads reprint.${mode}.tsv (CmdSaintExpress.java:59), got ${found[0][1].name}")
                }
                [aggregate_meta + [id: "saintexpress-${mode}".toString()], found ? found[0][1] : [], mode,
                 getToolArgs(configs, 'saintexpress'), requireToolField(configs, 'saintexpress', 'cmd_opts')]
            }
        }

    REPRINTTOSAINTINPUT(ch_for_saintexpress.map { meta, reprint_tsv, mode, _saint_cli, _cmd_opts -> [meta, reprint_tsv, mode] })

    SAINTEXPRESS(
        REPRINTTOSAINTINPUT.out.saint_input
            // no failOnMismatch: saint_input is optional, a mode without its reprint table has none (skipped)
            .join(ch_for_saintexpress.map { meta, _reprint_tsv, _mode, saint_cli, cmd_opts -> [meta, saint_cli, cmd_opts] })
    )

    // FragPipe's error for a skipped mode, also in the Nextflow log (skipped_mode.txt is a one-line task output)
    REPRINTTOSAINTINPUT.out.skipped.subscribe { item -> log.error(item[1].text.trim()) }  // item: [meta, skipped]
    ch_saint_skipped = REPRINTTOSAINTINPUT.out.skipped

    //
    // FPOP Quant (aggregate, CmdFpopQuant.java): after IonQuant / TMT-Integrator
    //
    ch_for_fpop = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'fpop') }
        .combine(ch_fpop_inputs)
        .map { configs, modpep, tmt_report, tmt_report_unmod ->
            [aggregate_meta, modpep, tmt_report, tmt_report_unmod, getToolArgs(configs, 'fpop')]
        }

    FPOP(ch_for_fpop)

    //
    // FPOP coADAPTr (aggregate, CmdFPOPcoadaptr.java): FPOP Modifications column in every
    // final psm.tsv; runs only when there are results (FragpipeRun.java:2163)
    //
    ch_coadaptr_tables = ch_psm_tables
        .mix(ch_diann_report.map { _meta, report -> ['diann', 'diann', [report.parent]] })

    // meta.fpop_source picks the publish path of the edited tables (conf/modules.config): the
    // path of the task whose unedited tables they replace
    ch_for_fpop_coadaptr = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'fpop_coadaptr') }
        .combine(ch_coadaptr_tables)
        .filter { configs, source, _layout, _dirs -> source == fpopCoadaptrSource(configs) }
        .map { configs, source, layout, dirs ->
            [aggregate_meta + [fpop_source: source], dirs, layout, getToolArgs(configs, 'fpop_coadaptr')]
        }

    FPOPCOADAPTR(ch_for_fpop_coadaptr)

    //
    // METAPROTEOMICS (aggregate, CmdMetaproteomics.java): one FP-Meta run over the whole
    // work dir, right after Percolator (FragpipeRun.java:2458)
    //
    ch_for_metaproteomics = ch_tool_configs
        .filter { configs -> shouldRunTool(configs, 'metaproteomics') }
        .combine(ch_meta_files.map { files -> [files] })
        .map { configs, files -> [aggregate_meta, files, getToolArgs(configs, 'metaproteomics')] }

    METAPROTEOMICS(
        ch_for_metaproteomics,
        ch_fasta,
        ch_taxon_files.map { names, _nodes -> names },
        ch_taxon_files.map { _names, nodes -> nodes },
        ch_decoy_tag
    )

    emit:
    saint_results    = SAINTEXPRESS.out.results_list           // channel: [ val(meta), path(list.txt) ]
    saint_skipped    = ch_saint_skipped                        // channel: [ val(meta), path(skipped_mode.txt) ] - modes FragPipe skips for a missing reprint table
    fpop_results     = FPOP.out.results_dir                    // channel: [ val(meta), path(results_dir) ]
    fpop_psms        = FPOPCOADAPTR.out.psms                  // channel: [ val(meta), path(psm.tsv) ]
    fpop_diann_report = FPOPCOADAPTR.out.diann_report         // channel: [ val(meta), path(report.tsv) ]
    taxonomy_results = METAPROTEOMICS.out.results_dir          // channel: [ val(meta), path(results_dir) ]
}
