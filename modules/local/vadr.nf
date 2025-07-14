process VADR {
  tag "$meta.id"
  label 'process_low'

  conda 'pkru22::vadr=1.6.3'
  container 'staphb/vadr:1.6.3'

  input:
  tuple val(meta), path(fasta)
  val(vadr_args)
  val(vadr_trim_args)

  output:
  tuple val(meta), path("${prefix}/*.vadr.pass.tbl"), optional: true, emit: feature_table
  tuple val(meta), path("${prefix}/*.vadr.pass.fa"), optional: true, emit: pass_fasta
  tuple val(meta), path("${prefix}/"), optional: true, emit: vadr_outdir
  path "versions.yml", optional: true, emit: versions

  script:
  def args = vadr_trim_args ?: ''
  def args2 = vadr_args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """

  if [ ! -s $fasta ]; then
        echo "Input FASTA file is empty. Skipping VADR process."
        exit 0
    fi
  
  fasta-trim-terminal-ambigs.pl ${args} $fasta > ${prefix}_trimmed.fasta

  if [ ! -s ${prefix}_trimmed.fasta ]; then
        echo "Trimmed FASTA is empty. Skipping VADR process."
        exit 0
    fi

  v-annotate.pl \\
    $args2 \\
    ${prefix}_trimmed.fasta \\
    ${prefix}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      vadr: \$(v-annotate.pl -h | perl -ne 'print "\$1\\n" if /^# VADR (\\d+\\.\\d+\\.\\d+)/')
  END_VERSIONS
  """
}


process VADR_SUMMARIZE_ISSUES {
  executor 'local'
  memory 100.MB

  input:
  path(vadr_output, stageAs: "input*/*")

  output:
  path('vadr-annotation-issues.txt'), emit: issues
  path('vadr-annotation-failed-sequences.txt'), emit: failed

  script:
  """
  cat input*/**/*.alt.list | awk 'NR == 1 || \$0 !~ /^#/' > vadr-annotation-issues.txt
  cat input*/**/*.fail.list > vadr-annotation-failed-sequences.txt
  """
}