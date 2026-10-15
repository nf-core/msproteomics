/*
 * Processes for tests/lib/read_only_inputs.nf.test: one write a task may make (a new file next to a staged
 * link, in the task dir) and the two it must not (next to the link's resolved target, inside a staged
 * directory). Under Fusion both land in the input bucket; tests/nextflow.config mounts inputs read-only
 * so they fail in the docker tests.
 */

process WRITE_NEXT_TO_LINK {
    container 'quay.io/nf-core/ubuntu:20.04'
    errorStrategy 'ignore'

    input:
    path input_file

    output:
    path "${input_file}.out", emit: out

    script:
    """
    cat ${input_file} > ${input_file}.out
    """
}

process WRITE_NEXT_TO_TARGET {
    container 'quay.io/nf-core/ubuntu:20.04'
    errorStrategy 'ignore'

    input:
    path input_file

    output:
    path 'done', emit: out

    script:
    """
    touch "\$(readlink -f ${input_file}).out"
    touch done
    """
}

process WRITE_INSIDE_STAGED_DIR {
    container 'quay.io/nf-core/ubuntu:20.04'
    errorStrategy 'ignore'

    input:
    path input_dir

    output:
    path 'done', emit: out

    script:
    """
    touch ${input_dir}/new_file
    touch done
    """
}

workflow READ_ONLY_INPUTS {
    take:
    input_file
    input_dir

    main:
    WRITE_NEXT_TO_LINK(input_file)
    WRITE_NEXT_TO_TARGET(input_file)
    WRITE_INSIDE_STAGED_DIR(input_dir)

    emit:
    next_to_link   = WRITE_NEXT_TO_LINK.out.out
    next_to_target = WRITE_NEXT_TO_TARGET.out.out
    inside_dir     = WRITE_INSIDE_STAGED_DIR.out.out
}
