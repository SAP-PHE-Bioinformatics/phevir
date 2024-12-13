
process NEXTCLADE {
  tag        "Clade Determination"
  label      "process_medium"
  container  'docker://nextstrain/nextclade:3.8.0'
  
  input:
  file(fasta)
  tuple val(dataset_name), path(dataset)

  output:
  path "nextclade/combined_${dataset_name}.csv", emit: nextclade_file
  path "nextclade/*", emit: results, optional: true
  path("nextclade/combined_${dataset_name}.tsv")           , emit: tsv

  tuple file("nextclade/${prefix}.aligned.fasta"), file("nextclade/${prefix}.nwk"), emit: prealigned, optional: true
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args ?: " "
  def files  = fasta.join(" ")
  def prefix = task.ext.prefix ?: "combined_${dataset_name}"
  """
    mkdir -p nextclade dataset logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    nextclade --version >> \$log
    nextclade_version=\$(nextclade --version)


    nextclade run ${args} \
      --input-dataset ${dataset} \
      --output-all=nextclade/ \
      --output-basename ${prefix} \\
      --jobs ${task.cpus} \
      $fasta \
      | tee -a \$log

    cp $fasta nextclade/${prefix}.fasta

    if [ -f "dataset/pathogen.json" ]
    then
      tag=\$(grep "tag" dataset/pathogen.json | grep tag | sed 's/\"//g' | sed 's/,//g' | awk '{print \$NF}')
    else
      tag="NA"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      nextclade: \$(nextclade --version | awk '{print \$NF}')
      tag: \$tag
      container: ${task.container}
    END_VERSIONS
  """
}