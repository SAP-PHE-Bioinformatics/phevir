process NEXTCLADE_DATASETGET {
  tag        "Downloading Dataset"
  label      "process_medium"
  container  'docker://nextstrain/nextclade:3.8.0'

  output:
  path "dataset", emit: dataset
  path "logs/${task.process}/${task.process}.${workflow.sessionId}.log"
  path "versions.yml", emit: versions

  shell:
  def args   = task.ext.args ?: "${params.nextclade_dataset}"
  """
    mkdir -p nextclade dataset logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    nextclade --version >> \$log
    nextclade_version=\$(nextclade --version)

    echo "Getting nextclade dataset for ${args}" | tee -a \$log
    nextclade dataset list | tee -a \$log

    nextclade dataset get --name ${args} --output-dir dataset

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      nextclade: \$(nextclade --version | awk '{print \$NF}')
      tag: \$(grep "tag" dataset/pathogen.json | grep tag | sed 's/\"//g' | sed 's/,//g' | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}