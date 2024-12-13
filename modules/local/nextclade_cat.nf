process NEXTCLADE_CAT {
    tag "$meta"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(cat_input)

    output:
    tuple val($baseName), path('NEXTCLADE_CLADE.tsv')     , emit: clades
    tuple val($baseName), path('NEXTCLADE_LINEAGE.tsv')   , optional: true, emit: lineage
    tuple val($baseName), path('*.nextclade_report.tsv')  , emit: nextclade_cat_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def baseName = cat_input.baseName
    def prefix = task.ext.prefix ?: "${baseName}"

    """
    python $projectDir/bin/nextclade_output_cat.py \\
        --id "${meta}.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}