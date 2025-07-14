process SUMMARY_IAV {
    tag "Create report"
    label 'process_medium'
    // using shiptv container since it has pandas, rich, typer installed
    conda '../../assets/test.yaml'

    container 'docker://iidjmay/python_odbc:latest'
    
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
    path("${runID}_IAV.csv"), emit: report
    path("${runID}_detailed_IAV.csv"), emit: detailed_report


    script:
    """

     summary_flu_pub.py --excel $blast \\
    --tsv $nextclade \\
    --run $runID \\
    --qc $qc --kraken2 $kraken2 --drug $drug_sensitivity_report --mut $mutation_report --typing $typing_report_tsv --irma $irma_consensus_qc


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
}