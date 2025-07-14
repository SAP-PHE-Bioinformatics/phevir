process DATAMASH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/datamash:1.1.0--0' :
        'biocontainers/datamash:1.1.0--0' }"
// docker pull quay.io/biocontainers/datamash:1.1.0--0

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $tsv | datamash \\
        min 3 mean 3 median 3 max 3 \\
        > ${prefix}_coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datamash: \$(echo \$(datamash --version 2>&1) | head -n 1 |cut -d " " -f 4)
    END_VERSIONS
    """
}