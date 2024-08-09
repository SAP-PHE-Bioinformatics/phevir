process IVAR_TRIM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/ivar:1.4.3--h43eeafb_0"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path('*.log'), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_trimmed"
    """
    ivar trim \\
        $args \\
        -i $bam \\
        -b $bed \\
        -p $prefix \\
        > ${prefix}.ivar.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*iVar version //; s/ .*\$//')
    END_VERSIONS
    """
}