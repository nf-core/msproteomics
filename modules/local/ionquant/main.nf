/*
 * IONQUANT: IonQuant quantification over FragPipe's per-group psm.tsv (CmdIonquant)
 *
 * ionquantRun builds filelist_ionquant.txt, the per-group TMT annotations and the option list as
 * CmdIonquant.configure (CmdIonquant.java:178-337); the module stages the group directories the way
 * FragPipe's work dir holds them, writes those files and runs IonQuant.
 */
include { fragpipeTools; requireFragpipeLicenseAgreement; writeTextFile; rawSymlinksToFiles; shellArgs; ionquantRun; ionquantStaging; ionquantFiles; ionquantCollect; fragpipeJarLines } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process IONQUANT {
    tag "${meta.id}"
    label 'process_high'

    input:
    // psm_dirs: PHILOSOPHERFILTER output directories, one per FragPipe group, each containing psm.tsv
    // spec_files: LC-MS files, staged in the task dir: FragPipe's --specdir is the LC-MS files' own directory
    //   (CmdIonquant.java:330-338), its work dir in a FragPipe run whose inputs sit there
    tuple val(meta), path(psm_dirs), path(spec_files)
    val(annotation_lines) // multi-plex TMT annotation lines (plex, channel, sample, sample_name, ...), [] when none
    val(config_cli)  // IonQuant options from the workflow parser, in CmdIonquant order (ionquantRun)
    val(modmasses)   // modmasses_ionquant.txt lines (list, may be empty), or false: no --modlist (O-Pair runs)
    val(root_group_dir) // psm_dirs entry holding FragPipe's blank group "" ('' when there is none)
    val(experiment_annotation) // experiment_annotation.tsv lines FragPipe writes for this pass ([] = none)
    path(glyco_masses)  // PTM-S ptm-shepherd-output/glyco_masses_list.txt for AppendFile ([] when AppendFile does not run)

    output:
    tuple val(meta), path("${prefix}")                    , emit: results_dir
    // FragPipe's blank group "" (its files in the work dir root), only when the run has one
    tuple val(meta), path("ionquant_root")                , emit: root_group, optional: true
    tuple val(meta), path("${prefix}/ion.tsv")            , emit: ions, optional: true
    tuple val(meta), path("${prefix}/combined_ion.tsv")   , emit: combined_ions, optional: true
    tuple val(meta), path("${prefix}/combined_protein.tsv"), emit: combined_protein, optional: true
    tuple val(meta), path("${prefix}/combined_peptide.tsv"), emit: combined_peptide, optional: true
    tuple val(meta), path("${prefix}/combined_modified_peptide.tsv"), emit: combined_modified_peptide, optional: true
    tuple val(meta), path("${prefix}/combined_site_*.tsv") , emit: combined_site, optional: true
    tuple val(meta), path("${prefix}/*.tsv")              , emit: all_tsv, optional: true
    tuple val(meta), path("experiment_annotation.tsv")    , emit: experiment_annotation, optional: true
    // IonQuant also writes the REPRINT tables of the work dir; SAINTexpress reads them (CmdSaintExpress.java:59)
    tuple val(meta), path("${prefix}/reprint.spc.tsv")    , emit: reprint_spc, optional: true
    tuple val(meta), path("${prefix}/reprint.int.tsv")    , emit: reprint_int, optional: true
    path "I_AGREE_FRAGPIPE_LICENSE_AGREEMENT", emit: license_agreement
    tuple val("${task.process}"), val('ionquant'), eval("cat .ionquant_version"), emit: versions_ionquant, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!task.memory) {
        error("${task.process} needs task.memory (JVM / Go memory limit)")
    }
    def mem = Math.max(1, task.memory.toGiga() - 2)
    def tools = fragpipeTools(task)
    def agreement = requireFragpipeLicenseAgreement(task)
    def groups = (psm_dirs instanceof List ? psm_dirs : [psm_dirs]).collect { d -> d.name }
    def run = ionquantRun(groups, root_group_dir, config_cli, annotation_lines, '{TASK_DIR}')
    def write_modlist = modmasses instanceof List
    // java -Djava.awt.headless=true -Xmx -Dlibs.bruker.dir -Dlibs.thermo.dir -cp <jfreechart>:<IonQuant>
    // ionquant.IonQuant --license <license> --threads (CmdIonquant.java:142-176): MSFragger's ext/ libraries,
    // the license FragPipe finds in tools/.. (FragpipeLocations.locateLicense)
    def argv = ['java', '-Duser.home={TASK_DIR}', '-Djava.awt.headless=true', "-Xmx${mem}G",
        "-Dlibs.bruker.dir=${tools.bruker_dir}", "-Dlibs.thermo.dir=${tools.thermo_dir}",
        '-cp', "${tools.jfreechart_jar}:${tools.ionquant_jar}", 'ionquant.IonQuant', '--license', tools.license,
        '--threads', "${task.cpus}"] + run.args + ['--filelist', '{TASK_DIR}/filelist_ionquant.txt'] +
        (write_modlist ? ['--modlist', '{TASK_DIR}/modmasses_ionquant.txt'] : [])
    """
    export HOME=\$(pwd)
    mkdir -p ${prefix}
    ${writeTextFile('I_AGREE_FRAGPIPE_LICENSE_AGREEMENT', agreement)}
    ${fragpipeJarLines(tools, ['ionquant', 'msfragger'])}
    echo "${tools.fragpipe_version}" > .fragpipe_version

${ionquantStaging(root_group_dir, run.named_groups)}
${ionquantFiles(run, write_modlist ? modmasses : null, experiment_annotation)}
    # AppendFile: PTM-S glycan masses appended to the work dir's modmasses_ionquant.txt, duplicates
    # (0.01 Da) skipped, with FragPipe's own AppendToFile (FragpipeRun.java:1822-1830, CmdAppendFile.java)
    ${glyco_masses ? "touch modmasses_ionquant.txt\n    java -Duser.home=\"\$PWD\" -cp \"${tools.fragpipe_lib}/*\" org.nesvilab.fragpipe.util.AppendToFile modmasses_ionquant.txt \"\$(pwd)/${glyco_masses}\"" : ''}

    ${rawSymlinksToFiles(spec_files)}

    CMD=(${shellArgs(argv)} ${args})
    printf '%s\\n' "\${CMD[*]}" > ${prefix}/ionquant.log
    "\${CMD[@]}" 2>&1 | tee -a ${prefix}/ionquant.log

${ionquantCollect(prefix, root_group_dir, run.named_groups)}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tools = fragpipeTools(task)
    def groups = (psm_dirs instanceof List ? psm_dirs : [psm_dirs]).collect { d -> d.name }
    def run = ionquantRun(groups, root_group_dir, config_cli, annotation_lines, '{TASK_DIR}')
    def write_modlist = modmasses instanceof List
    // The logged command omits the machine-dependent java, licence and --threads arguments
    def stub_argv = ['ionquant.IonQuant'] + run.args + ['--filelist', 'filelist_ionquant.txt'] + (write_modlist ? ['--modlist', 'modmasses_ionquant.txt'] : [])
    """
    mkdir -p ${prefix}
${rawSymlinksToFiles(spec_files)}
${ionquantStaging(root_group_dir, run.named_groups)}
${ionquantFiles(run, write_modlist ? modmasses : null, experiment_annotation)}
    CMD=(${shellArgs(stub_argv)} ${args})
    printf '%s\\n' "\${CMD[*]}" > ${prefix}/ionquant.log
    # IonQuant also writes the REPRINT tables of the work dir (reprint.spc.tsv, reprint.int.tsv in the golden FragPipe run)
    touch ion.tsv combined_ion.tsv combined_protein.tsv combined_peptide.tsv combined_modified_peptide.tsv reprint.spc.tsv reprint.int.tsv
${ionquantCollect(prefix, root_group_dir, run.named_groups)}
    touch I_AGREE_FRAGPIPE_LICENSE_AGREEMENT
    echo "${tools.fragpipe_version}" > .ionquant_version
    echo "${tools.fragpipe_version}" > .fragpipe_version
    """
}
