process MBG {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://fcyucn/fragpipe:24.0' :
        'docker.io/fcyucn/fragpipe:24.0' }"

    input:
    // MBG (Mass-Based Glycoproteomics) matching
    // psm_file: psm.tsv from PHILOSOPHER_FILTER
    // manifest_file: FragPipe manifest file
    // config_cli: MBG parameters as CLI string (coupled with sample inputs)
    tuple val(meta), path(psm_file), path(manifest_file), val(config_cli)
    path(residue_db)     // glycan residues database file (shared resource)
    path(glycan_mod_db)  // glycan modifications database file (shared resource)

    output:
    tuple val(meta), path("${prefix}")                          , emit: results_dir
    tuple val(meta), path("${prefix}/*_mbg*.tsv")               , emit: results, optional: true
    tuple val(meta), path("${prefix}/*_glycan*.tsv")            , emit: glycan_results, optional: true
    tuple val("${task.process}"), val('mbg'), eval("cat .mbg_version"), emit: versions_mbg, topic: versions
    tuple val("${task.process}"), val('fragpipe'), eval("cat .fragpipe_version"), emit: versions_fragpipe, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def ram = task.memory ? Math.max(1, task.memory.toGiga() - 2) : 30
    def tools_dir = task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'

    // MBG CLI parameters (from CmdMBGMatch.java):
    // --match: run matching mode
    // --psm: path to psm.tsv file
    // --manifest: path to manifest file
    // --toaddresiduals: glycan residues to add
    // --residuedb: path to glycan residues database
    // --glycanmoddb: path to glycan modifications database
    // --maxq: maximum q-value
    // --minpsms: minimum PSMs
    // --minglycans: minimum glycans
    // --fdr: FDR threshold
    // --mztol: m/z tolerance (ppm)
    // --rttol: retention time tolerance (min)
    // --imtol: ion mobility tolerance
    // --nopasef: disable PASEF (true/false)
    // --numthreads: number of threads
    // --runtmt: run TMT mode (true/false)
    // --expanddb: expand database (true/false)
    // --maxskips: max skips
    // --allowchimeric: allow chimeric spectra (true/false)

    """
    
    export HOME=\$(pwd)
    mkdir -p ${prefix}

    TOOLS_DIR="${tools_dir}"

    MBG_JAR=\$(find "\$TOOLS_DIR" -name 'MBG*.jar' -type f 2>/dev/null | head -1 || true)
    BATMASS_JAR=\$(find "\$TOOLS_DIR" -name 'batmass-io*.jar' -type f 2>/dev/null | head -1 || true)
    # IonQuant is optional for MBG
    IONQUANT_DEP_JAR=\$(find "\$TOOLS_DIR" -name 'IonQuant*.jar' -type f 2>/dev/null | head -1 || true)
    BRUKER_DIR="${tools_dir}/../ext/bruker"
    THERMO_DIR="${tools_dir}/../ext/thermo"
    _ver=\$(java -jar "\$MBG_JAR" --version 2>&1 | grep -oP 'Inference \\K[\\d.]+' | head -1 || true)

    echo "\$_ver" > .mbg_version
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)

    echo "\$_fp_ver" > .fragpipe_version
    CP="\$MBG_JAR:\$BATMASS_JAR"
    [ -n "\$IONQUANT_DEP_JAR" ] && CP="\$CP:\$IONQUANT_DEP_JAR"
    NATIVE_FLAGS=""
    [ -d "\$BRUKER_DIR" ] && NATIVE_FLAGS="\$NATIVE_FLAGS -Dlibs.bruker.dir=\$BRUKER_DIR"
    [ -d "\$THERMO_DIR" ] && NATIVE_FLAGS="\$NATIVE_FLAGS -Dlibs.thermo.dir=\$THERMO_DIR"

    # Run MBG matching
    export JAVA_OPTS="-Xmx${ram}G"
    CMD="java \$JAVA_OPTS \$NATIVE_FLAGS -cp \\"\$CP\\" com.mbg.MBG --match --psm ${psm_file} --manifest ${manifest_file} --residuedb ${residue_db} --glycanmoddb ${glycan_mod_db} --numthreads ${task.cpus} ${config_cli} ${args}"
    printf '%s\\n' "\$CMD" > ${prefix}/mbg.log
    eval "\$CMD" 2>&1 | tee -a ${prefix}/mbg.log

    # Move output files to results directory
    shopt -s nullglob
    for f in *_mbg*.tsv *_glycan*.tsv; do mv "\$f" ${prefix}/; done
    shopt -u nullglob

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/sample_mbg_results.tsv
    touch ${prefix}/sample_glycan_results.tsv
    touch ${prefix}/mbg.log
    TOOLS_DIR="${task.ext.fragpipe_tools_dir ?: '/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools'}"
    _fp_ver=\$(echo "\$TOOLS_DIR" | grep -oP 'fragpipe-\\K[0-9]+\\.[0-9]+[0-9.]*' | head -1 || true)
    echo "\$_fp_ver" > .mbg_version
    echo "\$_fp_ver" > .fragpipe_version

    """
}
