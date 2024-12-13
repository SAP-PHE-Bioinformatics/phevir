process RESISTANCE {
  tag "$fasta"
  label 'process_low'


  input:
  path(fasta)

  output:
  path("mutation_report.txt"), emit : mutation_report
  path("query_summary_report.csv"), emit : query_summary_report
  path("drug_sensitivity_report.tsv"), emit : drug_sensitivity_report

  script:
"""
  getresistance.sh $fasta 
"""
}
