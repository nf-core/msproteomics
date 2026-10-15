/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGPIPE_INFERENCE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Protein inference, FDR filtering and the multi-experiment report, as
    FragpipeRun.java:1670-1810 configures them.

    Modules:
    - PHILOSOPHERPROTEINPROPHET (all groups together, CmdProteinProphet)
    - PHILOSOPHERFILTER (DbAnnotate + filter + report per group, CmdPhilosopherFilter)
    - PHILOSOPHERIPROPHET (CmdIprophet) and PHILOSOPHERABACUS (CmdPhilosopherAbacus)

    Every tool sees FragPipe's work-dir layout: <group>/<file>, the blank group in the root.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PHILOSOPHERPROTEINPROPHET                                   } from '../../../modules/local/philosopherproteinprophet/main'
include { PHILOSOPHERFILTER } from '../../../modules/local/philosopherfilter/main'
include { PHILOSOPHERFILTER as PHILOSOPHERFILTER_LATER } from '../../../modules/local/philosopherfilter/main'
include { PHILOSOPHERIPROPHET                                         } from '../../../modules/local/philosopheriprophet/main'
include { PHILOSOPHERABACUS } from '../../../modules/local/philosopherabacus/main'
include { fpopCoadaptrMeta; shouldRunTool; getToolArgs; getToolField; getDecoyTag; requireToolField; fragpipeGroupDir; groupResultsDir; fragpipeLcmsOrderKey; compareFragpipeLcmsKeys; sortByRank; fragpipeGroupPath; fragpipeQuantLabels } from '../utils_nfcore_msproteomics_pipeline/fragpipe_utils'

workflow FRAGPIPE_INFERENCE {
    take:
    ch_pepxml        // channel: [ val(meta), path(pepxml(s)) ] - validated pepXMLs (sharedPepxmlFiles)
    ch_pepxml_group  // channel: [ val(meta), path(pepxml(s)) ] - every pepXML of a group dir (validated + PTMProphet)
    ch_fasta         // channel: path(fasta) - value channel for database
    ch_tool_configs  // channel: val(tool_configs_map) - parsed JSON config
    aggregate_meta   // val(meta) for aggregate outputs

    main:

    // Decoy tag (workflow database.decoy-tag) - single source for every decoy-aware tool
    ch_decoy_tag = ch_tool_configs.map { configs -> getDecoyTag(configs) }

    // sharedPepxmlFiles in FragPipe order (TreeMap of InputLcmsFile, InputLcmsFile.compareTo):
    // [ordered files, their work-dir locations, first group]
    ch_shared = ch_pepxml
        .toList()
        .map { items ->
            def ordered = items.sort(false) { a, b -> compareFragpipeLcmsKeys(fragpipeLcmsOrderKey(a[0]), fragpipeLcmsOrderKey(b[0])) }
            def entries = ordered.collectMany { meta, files -> sortByRank(files).collect { f -> [f, fragpipeGroupPath(meta.group, f.name)] } }
                .unique { entry -> entry[1] }
            [entries.collect { e -> e[0] }, entries.collect { e -> e[1] }, ordered ? ordered[0][0].group : null]
        }

    //
    // STEP 1: ProteinProphet over every group's pepXMLs (CmdProteinProphet.java:121-224)
    //
    ch_for_proteinprophet = ch_shared
        .combine(ch_tool_configs)
        .filter { files, _paths, _first, configs -> files && shouldRunTool(configs, 'proteinprophet') }
        .map { files, paths, _first, configs -> [aggregate_meta, files, paths, getToolArgs(configs, 'proteinprophet')] }

    PHILOSOPHERPROTEINPROPHET(ch_for_proteinprophet, ch_fasta)

    //
    // STEP 2: DbAnnotate + filter + report per group (FragpipeRun.java:1688-1780). Groups after the
    // first read the first group's workspace (--dbbin/--probin, CmdPhilosopherFilter.java:103-113).
    //
    ch_groups = ch_pepxml_group
        .toList()
        .map { items -> [items] }  // keep the list one element through combine
        .combine(ch_shared.map { _files, _paths, first -> first })
        .combine(ch_tool_configs)
        .filter { _items, first, configs -> first != null && shouldRunTool(configs, 'filter') }
        .flatMap { items, first, configs ->
            def by_group = items.groupBy { item -> item[0].group }
            def report_args = by_group.size() == 1 ? requireToolField(configs, 'filter', 'report_args_single_group') : getToolField(configs, 'filter', 'report_args')
            def dont_use = requireToolField(configs, 'filter', 'dont_use_protxml')
            // FreeQuant/Labelquant are followed by PHILOSOPHERREPORT, which writes the final tables
            def report_after_quant = shouldRunTool(configs, 'freequant') || shouldRunTool(configs, 'labelquant')
            by_group.collect { group, group_items ->
                def files = group_items.collectMany { item -> [item[1]].flatten() }.unique { f -> f.name }
                // FPOP coADAPTr replaces the group's final psm.tsv (filter's, or report's after quant);
                // FragPipe exports matched fragments into the final psm.tsv (FRAGPIPE_MATCHED_FRAGMENTS),
                // which then replaces the report's psm.tsv in the published results.
                def group_meta = [id: fragpipeGroupDir(group), group: group, report_after_quant: report_after_quant,
                                  export_matched_fragments: shouldRunTool(configs, 'exportmatchedfragments')]
                [fpopCoadaptrMeta(group_meta, configs, ['filter', 'report']), files, first,
                 getToolArgs(configs, 'filter'), dont_use, report_args]
            }
        }

    ch_protxml_for_filter = PHILOSOPHERPROTEINPROPHET.out.protxml.map { _meta, protxml -> [protxml] }
        .ifEmpty([[]])

    ch_filter_rows = ch_groups
        .combine(ch_protxml_for_filter)
        .map { meta, files, first, filter_args, dont_use, report_args, protxml ->
            def protxmls = [protxml].flatten().findAll { f -> f }
            if (!dont_use && !protxmls) {
                error("FRAGPIPE_INFERENCE: the filter uses the protXML but ProteinProphet produced none")
            }
            [meta, files, dont_use ? [] : protxmls[0], first, filter_args, dont_use, report_args]
        }
        .branch { row ->
            first: row[0].group == row[3]
            later: true
        }

    PHILOSOPHERFILTER(
        ch_filter_rows.first.map { meta, files, protxml, first, filter_args, dont_use, report_args ->
            [meta, files, protxml, [], first, filter_args, dont_use, report_args]
        },
        ch_fasta,
        ch_decoy_tag
    )

    PHILOSOPHERFILTER_LATER(
        ch_filter_rows.later
            .combine(PHILOSOPHERFILTER.out.workspace.map { _meta, workspace -> workspace })
            .map { meta, files, protxml, first, filter_args, dont_use, report_args, first_workspace ->
                [meta, files, protxml, first_workspace, first, filter_args, dont_use, report_args]
            },
        ch_fasta,
        ch_decoy_tag
    )

    ch_filter_files = PHILOSOPHERFILTER.out.results_files.mix(PHILOSOPHERFILTER_LATER.out.results_files)
    ch_workspaces = PHILOSOPHERFILTER.out.workspace.mix(PHILOSOPHERFILTER_LATER.out.workspace)

    //
    // STEP 3: iProphet and Abacus for more than one group without label-free IonQuant
    // (FragpipeRun.java:1782-1809).
    // Deliberate difference from FragPipe, documented in docs/usage.md (SAINTexpress):
    // Abacus --reprint also runs when SAINTexpress is on and no IonQuant run writes the REPRINT tables,
    // although FragPipe's own conditions would not run it; its outputs then exist only in the pipeline.
    //
    ch_report_groups = ch_workspaces
        .toList()
        .map { items -> [items] }  // keep the list one element through combine
        .combine(ch_tool_configs)
        .map { workspaces, configs ->
            // mapGroupsToProtxml is a TreeMap of LcmsFileGroup ordered by name (LcmsFileGroup.java:46-48)
            def ordered = workspaces.sort(false) { a, b -> a[0].group <=> b[0].group }
            def ionquant_reprint = shouldRunTool(configs, 'ionquant') || shouldRunTool(configs, 'ionquant_ms1') || shouldRunTool(configs, 'ionquant_isobaric')
            def abacus = ordered.size() > 1 && (shouldRunTool(configs, 'abacus') || (shouldRunTool(configs, 'saintexpress') && !ionquant_reprint))
            [ordered, abacus, abacus && shouldRunTool(configs, 'iprophet'), configs]
        }

    ch_for_iprophet = ch_report_groups
        .filter { _ordered, _abacus, iprophet, _configs -> iprophet }
        .combine(ch_shared)
        .map { _ordered, _abacus, _iprophet, configs, files, paths, _first -> [aggregate_meta, files, paths, getToolArgs(configs, 'iprophet')] }

    PHILOSOPHERIPROPHET(ch_for_iprophet)

    ch_for_abacus = ch_report_groups
        .filter { _ordered, abacus, _iprophet, _configs -> abacus }
        // Without ProteinProphet FragPipe still configures Abacus on combined.prot.xml of the work dir
        // (sharedMapGroupsToProtxml is refilled outside isRun(), FragpipeRun.java:1684-1685); Philosopher
        // then fails on the missing file when --protein reads it (philosopher lib/aba/pro.go).
        .combine(PHILOSOPHERPROTEINPROPHET.out.protxml.map { _meta, protxml -> [protxml] }.ifEmpty([[]]))
        .combine(PHILOSOPHERIPROPHET.out.pepxml.map { _meta, pepxml -> [pepxml] }.ifEmpty([[]]))
        .map { ordered, _abacus, _iprophet, configs, protxml_item, combined ->
            def protxml = [protxml_item].flatten().findAll { f -> f }
            def combined_pepxml = [combined].flatten().findAll { f -> f }
            def abacus_opts = requireToolField(configs, 'abacus', 'args').tokenize()
            if (abacus_opts.contains('--peptide') && !combined_pepxml) {
                error("FRAGPIPE_INFERENCE: abacus --peptide reads combined.pep.xml, which only iProphet writes (philosopher lib/aba/pep.go), and iProphet did not run")
            }
            if (abacus_opts.contains('--protein') && !protxml) {
                error("FRAGPIPE_INFERENCE: abacus --protein reads combined.prot.xml (philosopher lib/aba/pro.go: \"cannot find combined.prot.xml file\"), but ProteinProphet did not run")
            }
            def label = getToolField(configs, 'abacus', 'labels_channel_num')
            def plex = 0
            if (label) {
                if (fragpipeQuantLabels()[label] == null) {
                    error("FRAGPIPE_INFERENCE: unknown TMT label ${label} for abacus --plex")
                }
                plex = fragpipeQuantLabels()[label].size()
            }
            def abacus_args = requireToolField(configs, 'abacus', 'args')
            [aggregate_meta, ordered.collect { item -> item[1] }, ordered.collect { item -> item[0].group }, protxml ? protxml[0] : [],
             combined_pepxml ? combined_pepxml[0] : [], abacus_args, plex]
        }

    PHILOSOPHERABACUS(ch_for_abacus)

    ch_results_dir = ch_filter_files.map { meta, files -> [meta, groupResultsDir(files)] }

    emit:
    protxml          = PHILOSOPHERPROTEINPROPHET.out.protxml                  // channel: [ val(meta), path(combined.prot.xml) ]
    results_dir      = ch_results_dir                              // channel: [ val(meta), path(results_dir) ] per FragPipe group
    workspace        = ch_workspaces                               // channel: [ val(meta), path(.meta) ] per FragPipe group
    psms             = PHILOSOPHERFILTER.out.psms.mix(PHILOSOPHERFILTER_LATER.out.psms)          // channel: [ val(meta), path(psm.tsv) ]
    peptides         = PHILOSOPHERFILTER.out.peptides.mix(PHILOSOPHERFILTER_LATER.out.peptides)  // channel: [ val(meta), path(peptide.tsv) ]
    proteins         = PHILOSOPHERFILTER.out.proteins.mix(PHILOSOPHERFILTER_LATER.out.proteins)  // channel: [ val(meta), path(protein.tsv) ]
    ions             = PHILOSOPHERFILTER.out.ions.mix(PHILOSOPHERFILTER_LATER.out.ions)          // channel: [ val(meta), path(ion.tsv) ] optional
    pepxml           = ch_pepxml                                   // channel: [ val(meta), path(pepxml) ] - ProteinProphet doesn't modify pepXML
    combined_pepxml  = PHILOSOPHERIPROPHET.out.pepxml                         // channel: [ val(meta), path(combined.pep.xml) ]
    reprint_spc      = PHILOSOPHERABACUS.out.reprint_spc          // channel: [ val(meta), path(reprint.spc.tsv) ]
    reprint_int      = PHILOSOPHERABACUS.out.reprint_int          // channel: [ val(meta), path(reprint.int.tsv) ]
    combined_protein = PHILOSOPHERABACUS.out.combined_protein     // channel: [ val(meta), path(combined_protein.tsv) ]
    combined_peptide = PHILOSOPHERABACUS.out.combined_peptide     // channel: [ val(meta), path(combined_peptide.tsv) ]
    versions         = PHILOSOPHERPROTEINPROPHET.out.versions_philosopher
        .mix(PHILOSOPHERPROTEINPROPHET.out.versions_fragpipe)
        .mix(PHILOSOPHERFILTER.out.versions_philosopher)
        .mix(PHILOSOPHERFILTER.out.versions_fragpipe)
        .mix(PHILOSOPHERFILTER_LATER.out.versions_philosopher)
        .mix(PHILOSOPHERFILTER_LATER.out.versions_fragpipe)
        .mix(PHILOSOPHERIPROPHET.out.versions_philosopher)
        .mix(PHILOSOPHERIPROPHET.out.versions_fragpipe)
        .mix(PHILOSOPHERABACUS.out.versions_philosopher)
        .mix(PHILOSOPHERABACUS.out.versions_fragpipe)
}
