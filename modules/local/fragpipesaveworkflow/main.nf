/*
 * FRAGPIPESAVEWORKFLOW: Write FragPipe's saved UI state (fragpipe.workflow) for a run.
 *
 * FragPipe writes fragpipe.workflow into its work directory before any tool runs
 * (FragpipeRun.java:371-375). assets/fragpipe_ui_state/FragpipeUiState drives FragPipe's own
 * headless load and save code, so the file is FragPipe's: every UI key with its default, the
 * run values and the tool-version header. Needs the FragPipe image (--fragpipe_container).
 *
 * Input:
 *   - meta: Run metadata map
 *   - workflow_file: FragPipe .workflow file of the run, staged under input/: the output is fragpipe.workflow,
 *     and FragPipe would write through a staged link of that name into the user's file
 *   - manifest: fragpipe-files.fp-manifest of the run naming the staged spectra (spectra/<name>)
 *   - ui_state: directory with FragpipeUiState.class
 *   - database: the run database (database.db-path), staged under database/
 *   - workdir: absolute path FragPipe would run in (workdir)
 *   - spectra: the LC-MS files named in the manifest; FragPipe's headless manifest load stops when a
 *     listed file does not exist (TabWorkflow.showSkippedFiles), and staging them makes the container see them
 *
 * Output:
 *   - workflow: fragpipe.workflow
 *   - versions: Software versions (topic channel)
 */
include { fragpipeTools } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process FRAGPIPESAVEWORKFLOW {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(workflow_file, stageAs: 'input/*'), path(manifest), path(ui_state, stageAs: 'fragpipe_ui_state'), path(database, stageAs: 'database/*'), val(workdir), path(spectra, stageAs: 'spectra/*')

    output:
    tuple val(meta), path("fragpipe.workflow"), emit: workflow
    tuple val("${task.process}"), val('fragpipe'), eval("grep -m1 '^# FragPipe' fragpipe.workflow | sed 's/.* version //'"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def tools = fragpipeTools(task)
    """
    mkdir -p home tmp
    java -Djava.awt.headless=true -Duser.home=\$PWD/home -Djava.io.tmpdir=\$PWD/tmp ${args} \\
        -cp "${tools.fragpipe_lib}/*:fragpipe_ui_state" FragpipeUiState \\
        ${workflow_file} ${manifest} ${workdir} ${tools.tools_dir} ${database} fragpipe.workflow \\
        > fragpipe_ui_state.log 2>&1
    test -s fragpipe.workflow
    """

    stub:
    """
    printf '# Version info:\\n# FragPipe version 24.0\\n' > fragpipe.workflow
    """
}
