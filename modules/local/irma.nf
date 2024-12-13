process IRMA {
  tag "$meta.id"
  label 'process_high'

  // conda "bioconda::irma=1.0.2"
  // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
  //   container 'https://depot.galaxyproject.org/singularity/irma:1.0.2--pl5321hdfd78af_2'
  // } else {
  //   container 'quay.io/biocontainers/irma:1.0.2--pl5321hdfd78af_2'
  // }
  container 'docker://cdcgov/irma:v1.2.0'

  input:
  tuple val(meta), path(reads)
  val (irma_module)

  output:
  tuple val(meta), path("${meta.id}/"), emit: irma
  tuple val(meta), path("${meta.id}/*.bam")          , optional:true, emit: irma_bam
  tuple val(meta), path("${meta.id}/*.fasta")        , optional:true, emit: irma_fasta
  tuple val(meta), path("${meta.id}/*.vcf")          , optional:true, emit: irma_vcf
  tuple val(meta), path("*.irma.consensus.fasta")    , optional:true, emit: assembly
  tuple val(meta), path("*_LOW_ABUNDANCE.txt")       , optional:true, emit: failed_assembly
  tuple val(meta), path("*_HA.fasta")                , optional:true, emit: HA
  tuple val(meta), path("*_HA_FILE_NOT_FOUND.txt")   , optional:true, emit: failed_HA
  tuple val(meta), path("*_NA.fasta")                , optional:true, emit: NA
  tuple val(meta), path("*_NA_FILE_NOT_FOUND.txt")   , optional:true, emit: failed_NA
  tuple val(meta), path("*.irma_type.txt")           , emit: irma_type
  tuple val(meta), path("*.irma_subtype.txt")        , emit: irma_subtype
  tuple val(meta), path("*.irma.typing.tsv")         , emit: tsv
  tuple val(meta), path("${meta.id}.irma.consensus.fasta"), optional: true, emit: consensus
  tuple val(meta), path("${meta.id}.irma.majority_consensus.fasta"), optional: true, emit: majority_consensus
  path "*.irma.log", emit: log
  path "versions.yml", emit: versions

  script:
  def irma_config = "DEL_TYPE=\"NNN\"\nALIGN_PROG=\"BLAT\""
  def irma_log    = "${meta.id}.irma.log"
  def prefix = task.ext.prefix ?: "${meta.id}"
  def file = ''
  def found_coverage_files = ''
  """
  touch irma_config.sh
  echo 'SINGLE_LOCAL_PROC=${task.cpus}' >> irma_config.sh
  echo 'DOUBLE_LOCAL_PROC=${(task.cpus / 2).toInteger()}' >> irma_config.sh
  # default tmp in current working directory instead of defaulting to /tmp 
  # which may be restricted in size on HPC clusters
  echo 'ALLOW_TMP=1' >> irma_config.sh
  echo 'TMP=\$PWD' >> irma_config.sh
  if [ ${params.keep_ref_deletions} ]; then
    echo 'DEL_TYPE="NNN"' >> irma_config.sh
    echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
  fi

  IRMA $irma_module $reads $prefix

  if ls ${prefix}/amended_consensus/*.fa > /dev/null 2>&1; then
    cat ${prefix}/amended_consensus/*.fa > ${prefix}.irma.consensus.fasta
  else 
    echo "No consensus fasta due to a low abundance of read patterns per segment" > "${prefix}_LOW_ABUNDANCE"
    cat "${prefix}_LOW_ABUNDANCE" > "${prefix}_LOW_ABUNDANCE.txt"
  fi

  # Check and output IRMA type
  if [ -d "${prefix}" ] && [ -n "\$(ls -A "${prefix}"/*.fasta)" ]; then
      echo "Type_\$(basename \$(find "${prefix}" -name "*.fasta" | head -n1) | cut -d_ -f1)" > "${prefix}_IRMA_TYPE"
      cat "${prefix}_IRMA_TYPE" > "${prefix}.irma_type.txt"
  else
      echo "No IRMA type" > "${prefix}_IRMA_TYPE"
      cat "${prefix}_IRMA_TYPE" > "${prefix}.irma_type.txt"
  fi

  # Check for the presence of specific subtype files, process and consolidate as needed
  if [ -d "${prefix}" ] && [ -n "\$(ls -A ${prefix}/*HA_H*.fasta)" ]; then
      echo "\$(basename \$(find ${prefix} -name "*HA_H*.fasta" | head -n1 | rev | cut -d_ -f1 | rev))" > "${prefix}_HA_SUBTYPE"
  else
      echo "NoIRMAsubtype " > "${prefix}_HA_SUBTYPE"
  fi

  if [ -d "${prefix}" ] && [ -n "\$(ls -A ${prefix}/*NA_N*.fasta)" ]; then
      echo "\$(basename \$(find ${prefix} -name "*NA_N*.fasta" | head -n1 | rev | cut -d_ -f1 | rev))" > "${prefix}_NA_SUBTYPE"
  else
      echo "-NoIRMAsubtype" > "${prefix}_NA_SUBTYPE"
  fi

  if [ -s "${prefix}_HA_SUBTYPE" ] && [ -s "${prefix}_NA_SUBTYPE" ]; then
      cat "${prefix}_HA_SUBTYPE" "${prefix}_NA_SUBTYPE" > "${prefix}.subtype.txt"
      awk '{sub(".fasta","",\$1); printf \$1}' "${prefix}.subtype.txt" | sed 's/NoIRMAsubtype-NoIRMAsubtype/No IRMA subtype/' > "${prefix}.irma_subtype.txt"
  fi

  # Output IRMA typing tsv
  echo -e "Sample\tIRMA_type\tIRMA_subtype" > ${prefix}.irma.typing.tsv
  echo -e "${prefix}\t\$(cat ${prefix}.irma_type.txt)\t\$(cat ${prefix}.irma_subtype.txt)" >> ${prefix}.irma.typing.tsv

  if [ -f "${prefix}/amended_consensus/${prefix}_4.fa" ]; then
      cat "${prefix}/amended_consensus/${prefix}_4.fa" > "${prefix}_HA.fasta"
  else
      echo "No file found at ${prefix}/amended_consensus/${prefix}_4.fa" > "${prefix}_HA_FILE_NOT_FOUND"
      cat "${prefix}_HA_FILE_NOT_FOUND" > "${prefix}_HA_FILE_NOT_FOUND.txt"
  fi

  if [ -f "${prefix}/amended_consensus/${prefix}_6.fa" ]; then
      cat "${prefix}/amended_consensus/${prefix}_6.fa" > "${prefix}_NA.fasta"
  else
      echo "No file found at ${prefix}/amended_consensus/${prefix}_6.fa" > "${prefix}_NA_FILE_NOT_FOUND"
      cat "${prefix}_NA_FILE_NOT_FOUND" > "${prefix}_NA_FILE_NOT_FOUND.txt"
  fi

  if ls ${prefix}/tables/*-allAlleles.txt > /dev/null 2>&1; then
    irma-alleles2fasta -n "${prefix}" -i "${prefix}/tables" -o majority-consensus
    if ls majority-consensus/*.fasta > /dev/null 2>&1; then
      cat majority-consensus/*.fasta > ${prefix}.irma.majority_consensus.fasta
    fi
  fi

  ln -s .command.log $irma_log
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
     IRMA: \$(IRMA | head -n1 | sed -E 's/^Iter.*IRMA\\), v(\\S+) .*/\\1/')
  END_VERSIONS
  """
}
