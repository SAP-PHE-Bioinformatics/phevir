process mafft {
  tag           "Multiple Sequence Alignment"
  label         "process_high"
  //errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/mafft:7.520'

  when:
  task.ext.when == null || task.ext.when

  input:
  file(fasta)
  path(reference_genome)

  output:
  path "mafft/mafft_aligned.fasta", emit: msa, optional: true
  path "mafft/*", emit: files
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: " "
  def files  = fasta.join(" ")
  def prefix = task.ext.prefix ?: "mafft_aligned"
  """
    mkdir -p mafft logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    echo "mafft version:" >> \$log
    mafft --version 2>&1 >> \$log

    for fasta in ${files}
    do
      cat \$fasta >> mafft/ultimate.fasta
    done

    mafft --auto \
      ${args} \
      --thread ${task.cpus} \
      --add mafft/ultimate.fasta \
      ${reference_genome} \
      > mafft/${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      mafft: \$(mafft --version 2>&1 | awk '{print \$1 "_" \$2}' | head -n 1)
      container: ${task.container}
    END_VERSIONS
  """
}