process PANDEPTH {
    tag "$meta.id"
    label 'process_medium'
    container "docker://iidjmay/pandepth:2.25"

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("*.stat"), emit: stats
    path  "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    // Check if meta.segment exists and construct prefix accordingly
    def prefix = meta.segment ? "${meta.id}_${meta.segment}" : "${meta.id}"
    """
    pandepth \\
        -i ${bam_file} \\
        -o ${prefix} \\
        -t  ${task.cpus} \\
        -a \\
        ${args} \\
        ||  {
    # Generate a dummy file with zero values
    cat <<-EOF > ${prefix}.stat
	#Chr    Length  CoveredSite     TotalDepth      Coverage(%)     MeanDepth
	${meta.segment ?: 'unknown'}    0    0    0    0.00    0.00
	##RegionLength: 0    CoveredSite: 0       Coverage(%): 0.00     MeanDepth: 0.00
	EOF
    }
    if [ -f ${prefix}.chr.stat.gz ]; then
        gzip -d ${prefix}.chr.stat.gz
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandepth: \$(pandepth -v 2>&1 | head -n1 | sed 's/^.*pandepth //; s/ .*\$//')
    END_VERSIONS
    """
}