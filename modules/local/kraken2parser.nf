process KRAKENTOPMATCH {

    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container  'quay.io/biocontainers/pandas:1.5.2'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path(tsv), emit: top_match


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    kraken2_parser.py -i ${report} -o ${prefix}.tsv

    """
}