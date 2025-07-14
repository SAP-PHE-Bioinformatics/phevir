process SAMTOOLS_AMPLICONSTATS {
  tag        "${meta.id}"
  label      "process_single"
  container  'staphb/samtools:1.20'


  input:
  tuple val(meta), file(bam), file(index), val(bed_name), path(primer_bed)

  output:
  tuple val(meta.id), file("samtools_ampliconstats/${meta.id}_ampliconstats.txt"), emit: samtools_ampliconstats_files
  path "logs/${task.process}/${meta.id}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: " "
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
    mkdir -p samtools_ampliconstats logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    samtools ampliconstats ${args} \
      ${primer_bed} \
      ${bam} > samtools_ampliconstats/${prefix}_ampliconstats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}


process PLOT {
  tag           "${meta}"
  label         "process_single"
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  publishDir    path: params.outdir, mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/samtools:1.20'


  input:
  tuple val(meta), file(ampliconstats)

  output:
  path "samtools_plot_ampliconstats/${meta}*", emit: files
  path "logs/${task.process}/${meta}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args   ?: " " 
  def prefix = task.ext.prefix ?: "${meta}"
  """
    mkdir -p samtools_plot_ampliconstats/${prefix} logs/${task.process}
    log=logs/${task.process}/${prefix}.${workflow.sessionId}.log

    date > \$log
    samtools --version >> \$log

    plot-ampliconstats ${args} \
      samtools_plot_ampliconstats/${prefix} \
      ${ampliconstats}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      samtools: \$(samtools --version | head -n 1 | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}
