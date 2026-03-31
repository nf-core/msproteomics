process METAPROTEOMICS {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // Metaproteomics analysis for database optimization
    // project_dir: directory containing FragPipe results
    // config_cli: Metaproteomics parameters as CLI string (coupled with sample inputs)
    tuple val(meta), path(project_dir), val(config_cli)
    path(fasta)            // database FASTA file (shared resource)
    path(taxon_name_file)  // NCBI taxonomy names.dmp file (shared resource)
    path(taxon_node_file)  // NCBI taxonomy nodes.dmp file (shared resource)

    output:
    tuple val(meta), path("${prefix}")                              , emit: results_dir
    tuple val(meta), path("${prefix}/*_optimized.fasta")            , emit: optimized_fasta, optional: true
    tuple val(meta), path("${prefix}/*_taxonomy*.tsv")              , emit: taxonomy_results, optional: true
    tuple val(meta), path("${prefix}/*_metaproteomics*.tsv")        , emit: results, optional: true
    tuple val("${task.process}"), val('metaproteomics'), eval("cat .metaproteomics_version"), emit: versions_metaproteomics, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    // Metaproteomics CLI parameters (from CmdMetaproteomics.java):
    // DbOptimizer is the main command
    // --projectDir: path to project directory with FragPipe results
    // --fastaFile: path to FASTA database
    // --outdir: output directory
    // --decoyTag: decoy prefix (e.g., rev_)
    // --qvalue: q-value threshold
    // --deltaHyperscore: delta hyperscore threshold
    // --minPeptCntPerProt: minimum peptide count per protein
    // --minUniqPeptCntPerProt: minimum unique peptide count per protein
    // --minUniqPeptCnt: minimum unique peptide count
    // --taxonNameFile: NCBI taxonomy names.dmp file
    // --taxonNodeFile: NCBI taxonomy nodes.dmp file
    // --hostName: host organism name (quoted)
    // --iterations: number of iterations

    """

    export HOME=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    META_JAR=\$(find "\$TOOLS_DIR" -path '*/metaproteomics/FP-Meta*.jar' -type f 2>/dev/null | head -1 || true)
    _ver=\$(ls "\$TOOLS_DIR"/metaproteomics/FP-Meta-*.jar 2>/dev/null | grep -oP 'Meta-\\K[\\d.]+' | head -1 || true)

    echo "\$_ver" > .metaproteomics_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version

    # FP-Meta uses Java Files.walk() without FOLLOW_LINKS, so symlinked project
    # directories (as created by Nextflow staging) are not traversed.
    # Dereference symlinks by copying to a local directory.
    PROJECT_DIR_RESOLVED="project_resolved"
    cp -rL ${project_dir} "\$PROJECT_DIR_RESOLVED"

    export JAVA_OPTS="-Xmx${ram}G"
    CMD="java \${JAVA_OPTS:-} -jar \\"\$META_JAR\\" DbOptimizer --projectDir \$PROJECT_DIR_RESOLVED --fastaFile ${fasta} --outdir ${prefix} --taxonNameFile ${taxon_name_file} --taxonNodeFile ${taxon_node_file} ${config_cli} ${args}"
    printf '%s\\n' "\$CMD" > ${prefix}/metaproteomics.log
    eval "\$CMD" 2>&1 | tee -a ${prefix}/metaproteomics.log

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/sample_optimized.fasta
    touch ${prefix}/sample_taxonomy_results.tsv
    touch ${prefix}/sample_metaproteomics_results.tsv
    touch ${prefix}/metaproteomics.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .metaproteomics_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
