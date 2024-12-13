process SUMMARY {
    tag "Create report"
    label 'process_medium'
    // using shiptv container since it has pandas, rich, typer installed
    conda '../../assets/test.yaml'
    // // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    // //     container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.0--pyh5e36f6f_0'
    // // } else {
    // //     container 'quay.io/biocontainers/shiptv:0.4.0--pyh5e36f6f_0'
    // // }

    input:
    path(qc)
    path(irma_consensus_qc)
    path(typing_report_tsv)
    path(mutation_report) 
    path(drug_sensitivity_report)
    path(blast)
    path(nextclade)
    val(runID)
    path(kraken2)
   
    
    output:
    path("${runID}.csv"), emit: report

    script:
    """

     summary_flu.py --excel $blast \\
    --tsv $nextclade \\
    --run $runID \\
    --qc $qc --kraken2 $kraken2 --drug $drug_sensitivity_report --mut $mutation_report --typing $typing_report_tsv --irma $irma_consensus_qc


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
}