/*
 * FRAGPIPEPLAN: FragPipe's configured run for a workflow and manifest (no tool runs)
 *
 * Runs FragPipe headless with --write-plan --msfragger-per-file: FragPipe configures every tool as a real run
 * would, writes their config files and its run files (fragpipe.workflow, fragpipe-files.fp-manifest,
 * fragpipe.job) into its work dir fp_work/, logs every command (log_<time>.txt in fp_work/) and runs
 * none of them. Needs a FragPipe build with --write-plan and --msfragger-per-file (Nesvilab/FragPipe#2907;
 * no FragPipe release has them yet), passed with --fragpipe_container; with any other FragPipe the task
 * stops with that reason before running it.
 *
 * Input:
 *   - meta:            plan id map
 *   - raw_files:       LC-MS files, staged as links: FragPipe only checks that they exist, so unlike
 *                      FRAGPIPE_HEADLESS no .raw is copied (headlessStagingLines gets no raw files)
 *   - database:        FASTA (FragPipe counts entries and decoys)
 *   - workflow_text:   .workflow text (headlessWorkflowText sets the database path)
 *   - manifest_rows:   [file name, experiment, bioreplicate, data type, FragPipe group] per LC-MS file
 *   - annotation_rows: [group, channel, sample name] per TMT channel, [] for label-free
 *   - sources:         LC-MS file name -> original path (s3:// kept)
 *   - database_source: original database path
 *
 * Output:
 *   - workdir:           fp_work/, FragPipe's work dir (run log, configs, run files; PLANTOTOOLCONFIGS reads it)
 *   - manifest:          the manifest FragPipe read
 *   - portable_log:      the run log with the licence redacted, every input at its original path,
 *                        FragPipe's work dir @WORKDIR@ and its tools dir @TOOLS_DIR@ (bin/fragpipe_portable.py)
 *   - portable_workflow: the fragpipe.workflow FragPipe saves into its work dir, the same way
 *   - portable_manifest: the fragpipe-files.fp-manifest FragPipe saves into its work dir, the same way
 *   - versions_fragpipe: FragPipe version (versions topic)
 */
include { dryRunStubLogText; fragpipePlanCommand; fragpipeTools; headlessLayout; headlessStagingLines; planSourceMapText; shellArgs; writeTextFile } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process FRAGPIPEPLAN {
    tag "${meta.id}"
    label 'process_low'

    // No default container: the licensed FragPipe image comes from --fragpipe_container (process.container)

    input:
    tuple val(meta), path(raw_files, stageAs: "raw_files/*"), path(database), val(workflow_text), val(manifest_rows), val(annotation_rows), val(sources), val(database_source)

    output:
    tuple val(meta), path("fp_work", type: 'dir')               , emit: workdir
    tuple val(meta), path("manifest.fp-manifest")               , emit: manifest
    tuple val(meta), path("portable/log_*.txt")                 , emit: portable_log
    tuple val(meta), path("portable/fragpipe.workflow")        , emit: portable_workflow
    tuple val(meta), path("portable/fragpipe-files.fp-manifest"), emit: portable_manifest
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mem_gb = (task.memory.toGiga() * 0.9).intValue()
    def plan_run = fragpipePlanCommand(task, workflow_text, database.name, mem_gb)
    def tools = plan_run.tools
    def layout = headlessLayout(manifest_rows, annotation_rows)
    """
    # FragPipe's --help exits 1 (FragPipeMain.java:43-45); any other status is a failed java start. The plan
    # options must be listed exactly (a FragPipe release with --write-plan, Nesvilab/FragPipe#2907)
    set +e
    FRAGPIPE_HELP=\$(java -Duser.home=\$PWD -cp "${tools.fragpipe_lib}/*" org.nesvilab.fragpipe.FragPipeMain --help 2>&1)
    FRAGPIPE_HELP_STATUS=\$?
    set -e
    if [ \$FRAGPIPE_HELP_STATUS -ne 0 ] && [ \$FRAGPIPE_HELP_STATUS -ne 1 ]; then
        echo "FRAGPIPEPLAN: FragPipe --help failed (java exit \$FRAGPIPE_HELP_STATUS): \$FRAGPIPE_HELP" >&2
        exit 1
    fi
    for option in ${plan_run.options.join(' ')}; do
        if ! grep -qE -- "(^|[[:space:]])\$option([[:space:]]|\$)" <<< "\$FRAGPIPE_HELP"; then
            echo "FRAGPIPEPLAN: this FragPipe has no \$option (--fragpipe_planner plan needs a FragPipe release with --write-plan, Nesvilab/FragPipe#2907); pass such an image with --fragpipe_container" >&2
            exit 1
        fi
    done
    echo "${tools.fragpipe_version}" > .fragpipe_version

${headlessStagingLines([], layout)}
    ${writeTextFile('run.workflow', plan_run.workflow)}

    # user.home: the container user has no passwd entry (see FRAGPIPE_HEADLESS)
    JAVA_TOOL_OPTIONS="-Duser.home=\$PWD" ${shellArgs(plan_run.argv)} ${args}

    # Portable copies: staged inputs back to their original paths (exact map), FragPipe's work dir @WORKDIR@,
    # its tools dir @TOOLS_DIR@; stops if a task path is left
    ${writeTextFile('source_map.tsv', planSourceMapText(layout, manifest_rows, sources, database.name, database_source))}
    fragpipe_portable.py \\
        --task-dir \$PWD \\
        --tools-dir ${tools.tools_dir} \\
        --map source_map.tsv \\
        --outdir portable \\
        fp_work/log_*.txt fp_work/fragpipe.workflow fp_work/fragpipe-files.fp-manifest
    """

    stub:
    def layout = headlessLayout(manifest_rows, annotation_rows)
    """
${headlessStagingLines([], layout)}
    # a dry run with no commands: FragPipe's run log and run files in its work dir
    mkdir -p fp_work portable
    ${writeTextFile('fp_work/log_stub.txt', dryRunStubLogText(layout))}
    ${writeTextFile('fp_work/fragpipe.workflow', workflow_text)}
    ${writeTextFile('fp_work/fragpipe.job', "runtime-save\t{TASK_DIR}/fp_work/fragpipe-workflow.workflow\t{TASK_DIR}/fp_work/fragpipe-files.fp-manifest\t{TASK_DIR}/fp_work\t/stub/fragpipe/tools\t{TASK_DIR}/${database.name}\t1\t1\n")}
    cp manifest.fp-manifest fp_work/fragpipe-files.fp-manifest
    cp fp_work/log_stub.txt fp_work/fragpipe.workflow fp_work/fragpipe-files.fp-manifest portable/

    echo "${fragpipeTools(task).fragpipe_version}" > .fragpipe_version
    """
}
