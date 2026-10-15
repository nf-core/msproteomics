include { fragpipeJarLines; fragpipeTools } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { speclibgenCommand; taskDirArg; taskFileCommand; shellQuote; fragpipeToolVersions; javaXmxGb } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process SPECLIBGEN {
    tag "${meta.id}"
    label 'process_high'

    input:
    // FragPipe's SpecLibGen step (CmdSpecLibGen.java:114-216): one run of FragPipe's own
    // tools/speclib/gen_con_spec_lib.py in the group work dir, which holds, side by side,
    //   psm.tsv / peptide.tsv:     the group's Philosopher report (one group)
    //   pepxml_files:              the validated interact-*.pep.xml (or PTMProphet *.mod.pep.xml)
    //   spectra:                   the LC-MS files of filelist_speclibgen.txt (CmdSpecLibGen.java:172-202),
    //                              staged in the task dir (<work dir>/<name> in the file list)
    //   rt_cal_file / im_cal_file: user RT / IM calibration files, [] when the option is a keyword
    //   argv_parts:                gen_con_spec_lib.py argv values from PARSEFRAGPIPEWORKFLOW (speclibgen.argv_parts)
    tuple val(meta), path(psm, stageAs: 'psm.tsv'), path(peptide, stageAs: 'peptide.tsv'), path(pepxml_files), path(spectra), path(rt_cal_file, stageAs: 'rt_calibration/*'), path(im_cal_file, stageAs: 'im_calibration/*'), val(argv_parts)
    path(fasta)  // protein database (FragPipe's fastaPath, argv[1])

    output:
    tuple val(meta), path("library.tsv")                 , emit: library_tsv
    tuple val(meta), path("easypqp_files")               , emit: easypqp_files
    tuple val(meta), path("easypqp_rt_reference_run.tsv"), emit: rt_reference, optional: true
    // EasyPQP's ion mobility alignment, written for ion mobility data (diaPASEF) also with noIM: the reference run
    // and one alignment per run; gen_con_spec_lib.py deletes only the RT *.alignment_pkl files (:370-373)
    tuple val(meta), path("easypqp_im_reference_run.tsv"), emit: im_reference, optional: true
    tuple val(meta), path("easypqp_im_alignment_*.alignment_pkl"), emit: im_alignments, optional: true
    tuple val(meta), path("filelist_*.txt")              , emit: filelists
    tuple val(meta), path("speclibgen_argv.txt")         , emit: argv
    tuple val(meta), path("{irt.tsv,im.tsv,*_temp-psm.tsv,*.psmpkl,*.peakpkl,*_run_peaks.tsv,easypqp_lib_openswath.tsv}"), emit: intermediates, optional: true
    tuple val("${task.process}"), val('fragpipe-speclib'), val(fragpipeToolVersions().fragpipe_speclib), emit: versions_speclibgen, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (task.ext.args) {
        error("SPECLIBGEN: gen_con_spec_lib.py takes positional arguments only (CmdSpecLibGen.java:122-202); ext.args is not passed: ${task.ext.args}")
    }
    def ram = javaXmxGb(task, 2)
    def cmd = speclibgenCommand(argv_parts, fasta.name, spectra, rt_cal_file, im_cal_file, ram, task)
    """
    export HOME=\$(pwd)
    TASK_DIR=\$(pwd)
    ${fragpipeJarLines(fragpipeTools(task), ['msfragger'])}

    [ -f ${shellQuote(cmd.script)} ] || { echo ${shellQuote("ERROR: ${cmd.script} not found (the FragPipe tools directory set in nextflow.config)")} >&2; exit 1; }
    # binFragger.getBin() (CmdSpecLibGen.java:211)
    [ -f ${shellQuote(cmd.msfragger_jar)} ] || { echo ${shellQuote("ERROR: ${cmd.msfragger_jar} not found (the FragPipe tools directory set in nextflow.config)")} >&2; exit 1; }

    # FragPipe 24.0 pins fragpipe-speclib==${fragpipeToolVersions().fragpipe_speclib} (python/python_packages/requirements.txt)
    [ "\$(python3 -m pip show fragpipe-speclib | sed -n 's/^Version: //p')" = ${shellQuote(fragpipeToolVersions().fragpipe_speclib)} ] || { echo ${shellQuote("ERROR: fragpipe-speclib ${fragpipeToolVersions().fragpipe_speclib} is not installed (FragPipe 24.0 requirements.txt)")} >&2; exit 1; }
    echo ${shellQuote(cmd.fragpipe_version)} > .fragpipe_version

    # filelist_speclibgen.txt: one LC-MS path per line in LC-MS table order (CmdSpecLibGen.java:172-202)
    ${taskFileCommand('filelist_speclibgen.txt', cmd.filelist)}

    # Environment of the SpecLibGen process (CmdSpecLibGen.java:204-212)
    export PYTHONIOENCODING=utf-8
    export OMP_NUM_THREADS=${task.cpus}
    export FRAGPIPE_PARSE_MZML_CMD=${shellQuote(cmd.parse_mzml_cmd)}
    unset OPENMS_DATA_PATH

    # CmdSpecLibGen.java:123-202, cwd = group work dir
    SPECLIB_CMD=(${cmd.argv.collect { t -> taskDirArg(t) }.join(' ')})
    printf '%s\\n' "\${SPECLIB_CMD[@]}" > speclibgen_argv.txt
    "\${SPECLIB_CMD[@]}"
    """

    stub:
    def cmd = speclibgenCommand(argv_parts, fasta.name, spectra, rt_cal_file, im_cal_file, 30, task)
    """
    TASK_DIR=\$(pwd)
    ${taskFileCommand('filelist_speclibgen.txt', cmd.filelist)}
    SPECLIB_CMD=(${cmd.argv.collect { t -> taskDirArg(t) }.join(' ')})
    printf '%s\\n' "\${SPECLIB_CMD[@]}" > speclibgen_argv.txt

    touch library.tsv filelist_easypqp_library.txt easypqp_rt_reference_run.tsv
    # every output, as if the data had ion mobility (easypqp_im_*)
    touch easypqp_im_reference_run.tsv ${[spectra].flatten().collect { s -> shellQuote("easypqp_im_alignment_${s.baseName}.alignment_pkl") }.join(' ')}
    mkdir -p easypqp_files
    touch easypqp_files/cmds.txt
    echo ${shellQuote(cmd.fragpipe_version)} > .fragpipe_version
    """
}
