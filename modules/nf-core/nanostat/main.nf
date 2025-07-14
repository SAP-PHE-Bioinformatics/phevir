process NANOSTAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/nanoplot:1.41.6--pyhdfd78af_0' :
    //     'biocontainers/nanoplot:1.41.6--pyhdfd78af_0' }"
    container "quay.io/biocontainers/nanostat:0.1.5--py35_0"

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*.txt")                 , emit: txt
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_file = ("$ontfile".endsWith(".fastq.gz") || "$ontfile".endsWith(".fq.gz")) ? "--fastq ${ontfile}" :
        ("$ontfile".endsWith(".txt")) ? "--summary ${ontfile}" : ''
    """
    NanoStat \\
        $args \\
        -t $task.cpus \\
        $input_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanostat: \$(echo \$(NanoStat --version 2>&1) | sed 's/^.*NanoStat //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
  
    touch NanoStats.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
