process IRMA_SEGMENT_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    container 'docker://quay.io/gdsc/biopython-pandas-scipy:3.12.1'

    input:
    tuple val(meta), path(fasta_files)

    output:
    tuple val(meta), path("*_*.tsv") , emit: cov_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python $projectDir/bin/percent_cov_calc.py $fasta_files ${meta.id}
    """
}