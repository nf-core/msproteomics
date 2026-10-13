/*
 * EXPORTMATCHEDFRAGMENTS: FragPipe "Export matched fragments" (CmdExportMatchedFragments)
 *
 * FragPipe 24.0 runs, in its work dir (CmdExportMatchedFragments.java:58-70):
 *   java -cp tools/FP-PDV/FragAnnotator-1.0.0.jar:tools/batmass-io-<v>.jar FragmentsExportMainClass <workDir> <nThreads> r
 * The export reads fragpipe-files.fp-manifest, each group's psm.tsv (and protein.tsv header),
 * the spectra named in the manifest, and umpire-se.params (DIA-Umpire mode). It rewrites every
 * group's psm.tsv in place with three more columns: ions, ion_mz, ion_int.
 *
 * Input:
 *   - meta: aggregate metadata map
 *   - group_dirs: per-group Philosopher result directories (named by FragPipe group; the blank
 *     group "" uses the stand-in name root_group_dir)
 *   - spectra: spectra the manifest points to (manifest paths are ../spectra/<name>)
 *   - umpire_spectra: DIA-Umpire pseudo-spectra, staged in the FragPipe work dir like FragPipe
 *   - manifest: fragpipe-files.fp-manifest for this run
 *   - root_group_dir: group_dirs entry of the blank group ('' when there is none)
 *   - diaumpire, speclibgen: DIA-Umpire / SpecLibGen run (umpire-se.params; speclibgen.run-speclibgen
 *     in fragpipe.workflow, which makes the export read every group from the work-dir root)
 *   - export_psms: the psm.tsv files the export rewrites (matchedFragmentsExportPsms); the stub edits them
 *
 * Output:
 *   - results_dirs: results/<group dir>/ per group, psm.tsv with the matched fragments
 *   - psms: results/<group dir>/psm.tsv
 *   - log: console output of the export
 *   - versions: FragAnnotator version (from the jar name)
 */
include { fragpipeTools; matchedFragmentsStaging; matchedFragmentsCollect; fragpipeJarLines; shellArgs } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process EXPORTMATCHEDFRAGMENTS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(group_dirs, stageAs: 'groups/*'), path(spectra, stageAs: 'spectra/*'), path(umpire_spectra, stageAs: 'fragpipe/*')
    path(manifest, stageAs: 'fragpipe_manifest.tsv')
    val(root_group_dir)
    tuple val(diaumpire), val(speclibgen), val(export_psms)

    output:
    tuple val(meta), path("results/*", type: 'dir')      , emit: results_dirs
    tuple val(meta), path("results/*/psm.tsv")           , emit: psms
    tuple val(meta), path("export_matched_fragments.log"), emit: log
    tuple val("${task.process}"), val('fragannotator'), eval("cat .fragannotator_version"), emit: versions_fragannotator, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // FragPipe passes no options (CmdExportMatchedFragments.java); ext.args is appended after them
    def args = task.ext.args ?: ''
    def tools = fragpipeTools(task)
    def groups = (group_dirs instanceof List ? group_dirs : [group_dirs]).collect { d -> d.name.tokenize('/').last() }
    // TabRun.FRAG_ANNOTATOR and ToolingUtils.BATMASS_IO_JAR (batmass-io-1.36.6 of the release), in the tools folder.
    // FragAnnotator exits 0 when it finds no manifest or cannot read the results (ResultProcessor); FragPipe
    // then goes on with the psm.tsv files unchanged, and so does this module.
    """
    export HOME=\$(pwd)
    echo "1.0.0" > .fragannotator_version
    ${fragpipeJarLines(tools, ['batmass_io'])}

${matchedFragmentsStaging('fragpipe_manifest.tsv', groups, root_group_dir, diaumpire, speclibgen)}

    # FragPipe starts the export with the work dir as current directory (pb.directory(workDir))
    cd fragpipe
    java -Duser.home="\$PWD" -cp ${shellArgs(["${tools.fragannotator_jar}:${tools.batmass_io_jar}"])} FragmentsExportMainClass "\$PWD" ${task.cpus} r ${args} 2>&1 | tee ../export_matched_fragments.log
    cd ..

${matchedFragmentsCollect(groups, root_group_dir)}
    """

    stub:
    def groups = (group_dirs instanceof List ? group_dirs : [group_dirs]).collect { d -> d.name.tokenize('/').last() }
    """
${matchedFragmentsStaging('fragpipe_manifest.tsv', groups, root_group_dir, diaumpire, speclibgen)}

    # Stand-in for the export: the three matched-fragment columns appended to the psm.tsv files it rewrites
    for psm in ${export_psms.collect { f -> "'${f}'" }.join(' ')}; do
        [ -f "\$psm" ] || continue
        awk 'BEGIN { FS = OFS = "\\t" } NR == 1 { print \$0, "ions", "ion_mz", "ion_int"; next } { print \$0, "[]", "[]", "[]" }' "\$psm" > "\$psm.tmp"
        mv "\$psm.tmp" "\$psm"
    done
    echo "stub" > export_matched_fragments.log

${matchedFragmentsCollect(groups, root_group_dir)}
    echo "1.0.0" > .fragannotator_version
    """
}
