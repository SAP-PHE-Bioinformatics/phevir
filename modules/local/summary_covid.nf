process SUMMARY_COV {
    tag "${batch_id}-generate_summary"
    label 'process_medium'
    container 'docker://iidjmay/python_odbc:latest'

    input:
    path(QC) // Input coverage stats
    path(pangolin) // Pangolin results
    path(nextclade) // Nextclade results
    path(pangoCollapse) // Pango-collapse output
    path(fasta_files)
    val(run_id)

    output:
    path "${run_id}_covid_summary.csv"        // Output summary CSV
    path "${run_id}_passed_samples.txt"       // File with passed samples
    path "${run_id}_concatenated.fasta"       // Concatenated FASTA file
    path "${run_id}_gisaid.fasta"             // GISAID-formatted FASTA
    path "${run_id}_proforma.csv"             // Proforma file
    
    script:
    """

    mkdir -p fasta_dir
    mv $fasta_files fasta_dir/

    summary_cov_pub.py \
        --pandepth $QC \
        --pangolin $pangolin \
        --nextclade $nextclade \
        --pangocollapse $pangoCollapse \
        --run $run_id \
        --fasta_dir ./fasta_dir 
    """
}
