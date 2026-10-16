/*
 * FRAGPIPEHEADLESS: Run FragPipe in headless mode (all-in-one process)
 *
 * Runs the entire FragPipe pipeline (MSFragger, MSBooster, Percolator,
 * ProteinProphet, Philosopher Filter, IonQuant, etc.) in a single process
 * using FragPipe's built-in headless mode.
 *
 * This produces results identical to running FragPipe from the GUI.
 *
 * Input:
 *   - raw_files:            All raw data files (.d, .raw, .mzML)
 *   - database:             FASTA database file
 *   - workflow_file:        FragPipe .workflow configuration file
 *   - manifest_content:     Manifest lines file (filename\texperiment\tbioreplicate\tdata_type)
 *   - annotation_content:   TMT annotation file (group\tchannel\tsample_name), empty for LFQ
 *   - file_experiment_map:  LC-MS file name to FragPipe group file (filename\tgroup)
 *
 * Output:
 *   - all_results:       All FragPipe output files
 *   - combined_protein:  Combined protein report
 *   - combined_peptide:  Combined peptide report
 *   - combined_ion:      Combined ion report
 *   - manifest:          The manifest passed to FragPipe (--manifest)
 *   - annotations:       Per-group TMT annotation files (TMT only)
 *   - versions:          Software versions
 */
include { checkHeadlessRows; fragpipeTools; headlessLayout; headlessStagingLines; headlessWorkflowText; writeTextFile } from '../../../subworkflows/local/utils_nfcore_msproteomics_pipeline/fragpipe_utils_a'

process FRAGPIPEHEADLESS {
    tag "fragpipe_headless"
    label 'process_high'

    // No container directive: the FragPipe image given at runtime with --fragpipe_container (process.container)

    input:
    path raw_files, stageAs: "raw_files/*"
    path database
    val workflow_text    // the .workflow text FragPipe runs (headlessWorkflowText sets the database path)
    val manifest_rows    // [file name, experiment, bioreplicate, data type, FragPipe group] per LC-MS file
    val annotation_rows  // [group, channel, sample name] per TMT channel, [] for label-free

    output:
    path "results/**",                   emit: all_results
    path "results/combined_protein.tsv", emit: combined_protein, optional: true
    path "results/combined_peptide.tsv", emit: combined_peptide, optional: true
    path "results/combined_ion.tsv",     emit: combined_ion,     optional: true
    path "manifest.fp-manifest",         emit: manifest
    path "raw_files/**annotation.txt",   emit: annotations,      optional: true
    path "versions.yml",                 emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba (copied from modules/nf-core/diann/main.nf)
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FRAGPIPEHEADLESS module does not support Conda: it runs in the FragPipe image given with --fragpipe_container. Please use Docker / Singularity / Podman instead."
    }
    if (!task.memory) {
        error("${task.process}: task.memory is not set; FragPipe headless gets --ram and -Xmx from it")
    }
    def mem_gb = (task.memory.toGiga() * 0.9).intValue()
    def args = task.ext.args ?: ''
    def tools = fragpipeTools(task)
    checkHeadlessRows(raw_files, manifest_rows, annotation_rows)
    def layout = headlessLayout(manifest_rows, annotation_rows)
    """
    export JAVA_OPTS="-Xmx${mem_gb}G"
    mkdir -p results

${headlessStagingLines(raw_files, layout)}

    # FragPipe headless handles decoy generation internally via Philosopher.
    # Do NOT run philosopher here — it corrupts FASTA on Fusion-mounted filesystems.
    ${writeTextFile('run.workflow', headlessWorkflowText(workflow_text, database.name))}

    # Run FragPipe headless
    # FragPipe uses Gradle application plugin — launch via main class + classpath, not -jar.
    # user.home: the container user has no passwd entry, so the JVM resolves user.home to "?"
    # and writes ?/.cache, ?/.java into its cwd (results/ for the tools FragPipe spawns, which is
    # published). FragPipe builds its child java commands without user.home (cmd/Cmd*.java), so
    # only the environment reaches them: set it for this one command.
    JAVA_TOOL_OPTIONS="-Duser.home=\$PWD" java -Duser.home=\$PWD -Djava.awt.headless=true \${JAVA_OPTS} \\
        -cp "${tools.fragpipe_lib}/*" \\
        org.nesvilab.fragpipe.FragPipeMain \\
        --headless \\
        --workflow run.workflow \\
        --manifest manifest.fp-manifest \\
        --workdir \$(pwd)/results \\
        --threads ${task.cpus} \\
        --ram ${mem_gb} \\
        --config-tools-folder ${tools.tools_dir} ${args}

    # --help prints "FragPipePlus v24.0" on the first line and exits 1 (FragPipeMain.java:43-45): accept exactly that
    FRAGPIPE_HELP=\$(java -Duser.home=\$PWD -cp "${tools.fragpipe_lib}/*" org.nesvilab.fragpipe.FragPipeMain --help 2>&1) || test \$? -eq 1
    FRAGPIPE_VERSION=\$(printf '%s\\n' "\$FRAGPIPE_HELP" | grep -m1 -oP 'v\\K[0-9.]+')
    # printf, not a here-document: the staging lines start in column 0, so Nextflow strips no indentation
    # and an indented END marker would not end the here-document
    printf '"%s":\\n    fragpipe: "%s"\\n' "${task.process}" "\${FRAGPIPE_VERSION}" > versions.yml
    """

    stub:
    checkHeadlessRows(raw_files, manifest_rows, annotation_rows)
    def layout = headlessLayout(manifest_rows, annotation_rows)
    """
${headlessStagingLines(raw_files, layout)}
    mkdir -p results/sample1
    touch results/combined_protein.tsv
    touch results/combined_peptide.tsv
    touch results/combined_ion.tsv
    touch results/experiment_annotation.tsv
    ${annotation_rows ? "mkdir -p results/tmt-report && touch results/tmt-report/abundance_protein_None.tsv && printf 'tmtintegrator:\\n  log2transformed: true\\n' > results/tmt-integrator-conf.yml" : ''}
    touch results/sample1/psm.tsv
    touch results/sample1/protein.tsv

    printf '"%s":\\n    fragpipe: "24.0"\\n' "${task.process}" > versions.yml
    """
}
