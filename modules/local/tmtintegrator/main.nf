include { tmtIntegratorNames } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'
include { tmtIntegratorCommands; taskDirWord; taskFileCommand; shellQuote; fragpipeToolVersions } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_b'

process TMTINTEGRATOR {
    tag "${meta.id}"
    label 'process_high'

    input:
    // group_dirs: one directory per FragPipe group holding its psm.tsv/protein.tsv (IonQuant
    //             isobaric pass or Philosopher report), staged under input_results/
    // root_group_dir: the group_dirs entry holding FragPipe's blank group "" ('' when none)
    // config_text: tmt-integrator-conf.yml text from FragPipe's plan (output placeholder line)
    // msstats_annotation: MSstatsTMT_annotation.csv (msstatsTmtAnnotationLines, WRITETEXTFILE), [] when
    //                     Philosopher MSstats is off; copied into the work dir (CmdTmtIntegrator.java:105-140)
    tuple val(meta), path(group_dirs, stageAs: 'input_results/*'), val(root_group_dir), path(msstats_annotation, stageAs: 'msstats_input/*')
    val(config_text)

    output:
    tuple val(meta), path("${work_dir}/${report_dir}")                 , emit: results_dir
    tuple val(meta), path("${work_dir}/${report_dir}/abundance_*.tsv") , emit: abundance, optional: true
    tuple val(meta), path("${work_dir}/${report_dir}/ratio_*.tsv")     , emit: ratio, optional: true
    tuple val(meta), path("${work_dir}/${config_name}")                , emit: config
    tuple val(meta), path("${work_dir}/MSstatsTMT_annotation.csv")     , emit: msstats_annotation, optional: true
    tuple val("${task.process}"), val('tmtintegrator'), val(fragpipeToolVersions().tmtintegrator), emit: versions_tmtintegrator, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // FragPipe work dir: TMT-Integrator names the blank group's plex after the directory of its
    // psm.tsv, which FragPipe keeps in the work-dir root (CmdTmtIntegrator.java:294-296)
    def names = tmtIntegratorNames(task, prefix)
    work_dir = names.work_dir
    // FragPipe's second run for FPOP TMT writes tmt-report-unmod with tmt-integrator-conf-unmod.yml
    // (CmdTmtIntegrator.java:62, :149-152)
    report_dir = names.report_dir
    config_name = names.config_name
    if (!task.memory) {
        error("${task.process} needs task.memory (JVM / Go memory limit)")
    }
    def mem = Math.max(1, task.memory.toGiga() - 2)
    def cmds = tmtIntegratorCommands(group_dirs, root_group_dir, config_text, work_dir, report_dir, task)
    def msstats = msstats_annotation ? "cp msstats_input/MSstatsTMT_annotation.csv ${shellQuote(work_dir)}/" : ''
    """
    export HOME=\$(pwd)
    TASK_DIR=\$(pwd)
    TMTINTEGRATOR_JAR=${shellQuote(cmds.jar)}
    [ -f "\$TMTINTEGRATOR_JAR" ] || { echo "ERROR: \$TMTINTEGRATOR_JAR not found (the FragPipe tools directory set in nextflow.config)" >&2; exit 1; }
    echo ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version

    # FragPipe's work dir: the blank group at its top, named groups below it (CmdTmtIntegrator.java:144-170, 288-296)
    ${cmds.stage.join('\n    ')}
    ${taskFileCommand("${work_dir}/${config_name}", cmds.config)}
    ${msstats}

    # Config check, then java -Xmx -jar <TMT-Integrator> <conf> <psm.tsv per group>, run in the
    # work dir (CmdTmtIntegrator.java:254-300)
    cd ${work_dir}
    java -Duser.home="\$PWD" -jar "\$TMTINTEGRATOR_JAR" "\$TASK_DIR"/${shellQuote("${work_dir}/${config_name}")} --ValParam
    CMD=(java -Duser.home="\$PWD" "-Xmx${mem}G" -jar "\$TMTINTEGRATOR_JAR" "\$TASK_DIR"/${shellQuote("${work_dir}/${config_name}")} ${cmds.psm_files.collect { p -> taskDirWord(p) }.join(' ')} ${args})
    printf '%s\\n' "\${CMD[*]}" > ../tmt-integrator.log
    "\${CMD[@]}" 2>&1 | tee -a ../tmt-integrator.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def names = tmtIntegratorNames(task, prefix)
    work_dir = names.work_dir
    report_dir = names.report_dir
    config_name = names.config_name
    def cmds = tmtIntegratorCommands(group_dirs, root_group_dir, config_text, work_dir, report_dir, task)
    def msstats = msstats_annotation ? "cp msstats_input/MSstatsTMT_annotation.csv ${shellQuote(work_dir)}/" : ''
    """
    TASK_DIR=\$(pwd)
    ${cmds.stage.join('\n    ')}
    ${taskFileCommand("${work_dir}/${config_name}", cmds.config)}
    ${msstats}
    printf '%s\\n' "TMT-Integrator ${config_name} "${cmds.psm_files.collect { p -> taskDirWord(p) }.join('" "')} > tmt-integrator.log
    touch ${work_dir}/${report_dir}/abundance_gene_None.tsv
    touch ${work_dir}/${report_dir}/abundance_protein_None.tsv
    touch ${work_dir}/${report_dir}/ratio_gene_None.tsv
    touch ${work_dir}/${report_dir}/ratio_protein_None.tsv
    echo ${shellQuote(cmds.fragpipe_version)} > .fragpipe_version
    """
}
