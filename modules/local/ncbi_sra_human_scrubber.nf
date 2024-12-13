process NCBI_SRA_HUMAN_SCRUBBER {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-human-scrubber:2.0.0--hdfd78af_0':
        'quay.io/hdc-workflows/sra-human-scrubber:2.0.0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_dehosted.fastq.gz"), emit: reads
    tuple val(meta), path('*.SPOTS_REMOVED.txt')       , emit: spots_removed
    path 'versions.yml'                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    

    """
    # Unzip the first read file if it's gzipped
    if [[ "$reads" == *.gz ]]; then
        echo "Decompressing $reads"
        gunzip -c $reads > reads.fastq
    fi


    # Run the scrubbing tool on each read and capture the count of masked reads
    scrub.sh reads.fastq |& tail -n1 | awk -F" " '{print \$1}' > SPOTS_REMOVED
    
    # Compress the dehosted/cleaned reads into gzipped fastq format
    gzip reads.fastq.clean -c > ${meta.id}_dehosted.fastq.gz


    if [ -f "SPOTS_REMOVED" ]; then
        cat "SPOTS_REMOVED" > "${meta.id}.SPOTS_REMOVED.txt"
    fi


    ## Calculate the total spots masked across both forward and reverse reads
    if [ -f "${meta.id}.SPOTS_REMOVED.txt" ]; then
        awk '{s+=\$1} END {print s}' ${meta.id}.SPOTS_REMOVED.txt > TOTAL_SPOTS_REMOVED
    fi
   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi_sra_human_scrubber: 2.0.0
    END_VERSIONS
    """
}