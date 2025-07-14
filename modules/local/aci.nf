process ACI{
    tag        "${meta.id}"
    label      "process_high"
    container  'staphb/aci:1.4.20240116'
    errorStrategy { task.attempt < 3 ? 'retry' : 'ignore'}

    input:
    tuple val(meta), file(bam), file(bai)
    path(bed)

    output:
    path "aci/${meta.id}/${meta.id}_amplicon_depth.csv", emit: cov, optional: true
    path "aci/${meta.id}/${meta.id}_amplicon_depth.png", emit: for_multiqc, optional: true
    path "aci/${meta.id}/*", emit: everything
    path "logs/${task.process}/*.log"
    path "versions.yml", emit: versions
  
    shell:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
        mkdir -p aci/${prefix} logs/${task.process}
        log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

        # time stamp + capturing tool versions
        date > \$log
        aci --version >> \$log

        aci ${args} \
            --bam ${bam} \
            --bed ${bed} \
            --threads ${task.cpus} \
            --out aci/${prefix} \
            | tee -a \$log
        
        if [ -f "aci/${prefix}/amplicon_depth.csv" ] ; then cp aci/${prefix}/amplicon_depth.csv aci/${prefix}/${prefix}_amplicon_depth.csv ; fi
        if [ -f "aci/${prefix}/amplicon_depth.png" ] ; then cp aci/${prefix}/amplicon_depth.png aci/${prefix}/${prefix}_amplicon_depth.png ; fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            aci: \$(aci --version | awk '{print \$NF}')
            container: ${task.container}
        END_VERSIONS
    """
}