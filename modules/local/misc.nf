process CAT_DB {
    tag "$fasta1 - $fasta2"

    executor 'local'
    memory 100.MB

    input:
    path(fasta1)
    path(fasta2)

    output:
    path("influenza_db.fasta"), emit: fasta

    script:
    """
    cp $fasta1 influenza_db.fasta
    echo >> influenza_db.fasta
    cat $fasta2 >> influenza_db.fasta
    """
}

process CAT_CONSENSUS {
  tag "$sample"
  conda 'bioconda::shiptv=0.4.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.0--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.0--pyh5e36f6f_0'
  }

  input:
  tuple val(sample), path(consensus)

  output:
  tuple val(sample), path('*.consensus.blastn.fasta'), emit: fasta
  path('*.consensus.fasta'), emit: consensus_fasta
  path "versions.yml" , emit: versions

  script:
  """
  cat_consensus_sequences.py \\
    --sample-name $sample \\
    --output1-fasta ${sample}.consensus.fasta \\
    --output2-fasta ${sample}.consensus.blastn.fasta \\
    $consensus

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     python: \$(python --version | sed 's/Python //g')
  END_VERSIONS
  """
}

// pass in a ch of samples from a .collect and concatenate them into one fasta output file
//input is a list of tuples with sample name and fasta file, runID used to call new file
// output is the fasta concatenated file with all the samples
process CAT_FASTA {
  tag "$runID"
  conda 'bioconda::shiptv=0.4.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/shiptv:0.4.0--pyh5e36f6f_0'
  } else {
    container 'quay.io/biocontainers/shiptv:0.4.0--pyh5e36f6f_0'
  }

  input:
  val(runID)
  path(dir)

  output:
  tuple val(runID), path('*.fasta'), emit: fasta
  path "versions.yml" , emit: versions

  script:
  """
  cat $dir/*.fasta >> ${runID}.fasta

  END_VERSIONS
  """
}

def fluPrefix(sample, segment, ref_id) {
    return "${sample}.Segment_${segment}.${ref_id}"
}
