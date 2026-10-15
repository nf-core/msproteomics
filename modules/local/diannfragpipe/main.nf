/*
 * DIANNFRAGPIPE: FragPipe's DIA-NN quantification step and its post-processing (CmdDiann)
 *
 * Input:
 *   - meta: Run metadata map
 *   - library: SpecLibGen library.tsv or the diann.library file
 *   - spectra: DIA LC-MS files
 *   - psm / peptide / protein: FragPipe reports DIA-NN's post-processing reads
 *   - manifest: fragpipe-files.fp-manifest
 *   - diann_config: DIA-NN settings (JSON) from the workflow parser
 *   - fasta: protein database
 *
 * Output:
 *   - results_files / report / pr_matrix / pg_matrix / speclib / msstats: dia-quant-output/ files
 *   - filelist / argv: filelist_diann.txt and the DIA-NN argv
 *   - versions_diann / versions_fragpipe: versions (topic: versions)
 */
include { fragpipeTools; diannArgv; diannFilelistText; diannInputDirLines; writeTextFile; shellArgs; fragpipeJarLines; rawLinkCopies } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { javaXmxGb } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

// Long work dirs: DIA-NN 1.8.2 Beta 8 aborts on a run path of 256+ characters, so the runs are then read through a
// short /tmp link (diannInputDirLines) and the outputs record /tmp/diann.XXXXXX/s/<run> as the runs' location.
process DIANNFRAGPIPE {
    tag "${meta.id}"
    label 'process_high'

    input:
    // FragPipe's DIA-Quant step (CmdDiann.configure, CmdDiann.java:100-621) in the group work dir:
    //   library:               SpecLibGen's library.tsv, or the workflow's diann.library file
    //   spectra:               the DIA / DIA-Quant LC-MS files (diannSpectraFiles), largest first (diannSizeOrder),
    //                          staged in the task dir: filelist_diann.txt and the report run names then read
    //                          <work dir>/<run>.d, as FragPipe's with the LC-MS files beside its work dir
    //   psm, peptide, protein: the group's Philosopher reports; Propagation and DiannToMsstats read psm.tsv,
    //                          FragReporter reads peptide.tsv and protein.tsv next to it (without them it
    //                          falls back to "PSM data only" and its abundance_* reports differ)
    //   manifest:              fragpipe-files.fp-manifest (DiannToMsstats, FragReporter)
    //   diann_config:          JSON of the parsed diann entry (args, args_post and run settings)
    tuple val(meta), path(library, stageAs: 'library_input/*'), path(spectra), path(psm, stageAs: 'psm.tsv'), path(peptide, stageAs: 'peptide.tsv'), path(protein, stageAs: 'protein.tsv'), path(manifest, stageAs: 'fragpipe-files.fp-manifest'), val(diann_config)
    path(fasta)  // FragPipe's fastaFile (--relaxed-prot-inf --fasta)

    output:
    tuple val(meta), path("dia-quant-output/*")                            , emit: results_files
    tuple val(meta), path("dia-quant-output/report.tsv")                   , emit: report
    tuple val(meta), path("dia-quant-output/report.pr_matrix.tsv")         , emit: pr_matrix
    tuple val(meta), path("dia-quant-output/report.pg_matrix.tsv")         , emit: pg_matrix
    tuple val(meta), path("dia-quant-output/report-tsv.speclib")           , emit: speclib, optional: true
    tuple val(meta), path("dia-quant-output/msstats.csv")                  , emit: msstats, optional: true
    tuple val(meta), path("filelist_diann.txt")                            , emit: filelist
    tuple val(meta), path("diann_argv.txt")                                , emit: argv
    tuple val("${task.process}"), val('diann'), eval("cat .diann_version"), emit: versions_diann, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba (copied from modules/nf-core/diann/main.nf)
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DIANNFRAGPIPE module does not support Conda: it runs in the FragPipe image given with --fragpipe_container. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def cfg = new groovy.json.JsonSlurper().parseText(diann_config)
    def ram = javaXmxGb(task, 2)
    def tools = fragpipeTools(task)
    if (!cfg.library && library.toString().tokenize('/').last() != 'library.tsv') {
        error("DIANNFRAGPIPE ${meta.id}: without diann.library DIA-NN reads SpecLibGen's library.tsv (CmdDiann.java:313-326), got ${library}")
    }
    def inputs = [spectra].flatten()
    def argv = diannArgv(cfg, tools.diann, task.cpus, library, fasta)
    def plex_labels = [['--light', cfg.light], ['--medium', cfg.medium], ['--heavy', cfg.heavy]].findAll { pair -> pair[1] }.collectMany { pair -> pair }
    // DIA-NN is FragPipe 24.0's bundled binary (Diann.DIANN_LINUX, Diann.java:56-58). Diann.getLDPRELOAD
    // (Diann.java:161-190) sets no LD_PRELOAD when the binary runs bare, as it does in the FragPipe image.
    // Diann.validate ignores the exit status of "--version" (Diann.java:117-131); the banner check decides.
    """
    ${rawLinkCopies(spectra)}
    export HOME=\$(pwd)
    WD=\$(pwd)
    echo "${tools.fragpipe_version}" > .fragpipe_version

    # The parser emits DIA-NN 1.8.2 Beta 8 arguments (NoteConfigDiann.compareVersion branches)
    "${tools.diann}" --version > diann_version.log 2>&1 || echo "exit status \$?" >> diann_version.log
    grep -oiP 'DIA-NN \\K[\\w\\s.-]+(?= \\(Data-Independent Acquisition by Neural Networks\\))' diann_version.log | head -1 > .diann_version
    if ! grep -qix '1.8.2 Beta 8' .diann_version; then
        echo "ERROR: ${tools.diann} reports version '\$(cat .diann_version)', FragPipe 24.0 runs DIA-NN 1.8.2 Beta 8" >&2
        exit 1
    fi

    ${fragpipeJarLines(tools, ['batmass_io'])}
    ${cfg.library ? '' : 'ln -s library_input/library.tsv library.tsv'}
    mkdir -p dia-quant-output

    # plexDIA library (CmdDiann.java:142-197)
    ${cfg.run_plex ? "java -Duser.home=\"\$PWD\" -cp ${shellArgs(["${tools.batmass_io_jar}:${tools.fragpipe_lib}/*"])} org.nesvilab.fragpipe.tools.diann.PlexDiaHelper --threads ${task.cpus} ${shellArgs(plex_labels)} --library library.tsv --out library_2.tsv" : ''}

    ${diannInputDirLines(inputs)}
    ${writeTextFile('filelist_diann.txt', diannFilelistText(inputs))}

    DIANN_CMD=(${shellArgs(argv.head)} ${args} ${shellArgs(argv.tail)})
    printf '%s\\n' "\${DIANN_CMD[@]}" > diann_argv.txt
    "\${DIANN_CMD[@]}"

    # Delete the .quant files next to the inputs (CmdDiann.java:411-417)
    rm -f ${shellArgs(inputs.collect { f -> "{TASK_DIR}/${f}.quant" })}

    # Move DIA-NN's .speclib to dia-quant-output/report-tsv.speclib (CmdDiann.java:436-438)
    java -Duser.home="\$PWD" -cp "${tools.fragpipe_jar}" org.nesvilab.utils.FileMoveWithExtension "\$WD" .speclib "\$WD/dia-quant-output/report-tsv.speclib"

    cd dia-quant-output
    # Propagate protein/gene information (CmdDiann.java:460-488)
    ${cfg.run_plex ? '' : "java -Duser.home=\"\$PWD\" -cp ${shellArgs(["${tools.batmass_io_jar}:${tools.fragpipe_lib}/*"])} org.nesvilab.fragpipe.tools.diann.Propagation \"\$WD\" ${shellArgs([cfg.light ?: '-', cfg.medium ?: '-', cfg.heavy ?: '-'])}"}
    # MSstats input (CmdDiann.java:490-517)
    ${!cfg.run_plex && cfg.generate_msstats ? "java -Duser.home=\"\$PWD\" -cp \"${tools.fragpipe_lib}/*\" org.nesvilab.fragpipe.tools.diann.DiannToMsstats report.tsv ./ \"\$WD/psm.tsv\" ${cfg.qvalue} ${cfg.run_specific_protein_qvalue ? cfg.qvalue : '1'} ${cfg.qvalue} ${cfg.qvalue} \"\$WD/fragpipe-files.fp-manifest\"" : ''}
    # FragReporter site / peptide reports (CmdDiann.java:519-547)
    java -Duser.home="\$PWD" -Xmx${ram}G -jar "${tools.fragreporter_jar}" --pr "\$WD/dia-quant-output/report.tsv" --exp-ann "\$WD/fragpipe-files.fp-manifest" \\
        --out-dir "\$WD/dia-quant-output" --mod-tag ${shellArgs([cfg.mod_tag])} --min-site-prob ${cfg.min_site_prob} --level ${shellArgs([cfg.report_levels])}
    cd "\$WD"
    """

    stub:
    def args = task.ext.args ?: ''
    def cfg = new groovy.json.JsonSlurper().parseText(diann_config)
    def tools = fragpipeTools(task)
    def inputs = [spectra].flatten()
    def argv = diannArgv(cfg, tools.diann, task.cpus, library, fasta)
    """
    ${diannInputDirLines(inputs)}
    ${writeTextFile('filelist_diann.txt', diannFilelistText(inputs))}
    DIANN_CMD=(${shellArgs(argv.head)} ${args} ${shellArgs(argv.tail)})
    printf '%s\\n' "\${DIANN_CMD[@]}" > diann_argv.txt
    mkdir -p dia-quant-output
    touch dia-quant-output/report.tsv dia-quant-output/report.pr_matrix.tsv dia-quant-output/report.pg_matrix.tsv dia-quant-output/report-tsv.speclib
    ${cfg.generate_msstats ? 'touch dia-quant-output/msstats.csv' : ''}
    echo "1.8.2 Beta 8" > .diann_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}
