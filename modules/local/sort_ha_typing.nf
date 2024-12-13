
process SORT_H_TYPING {
    memory {// Dynamically determine how much memory is required for this task based on 
    // overall size of tabular blastn inputs. For a single input, allocate 2GB
    if (blastn_results instanceof nextflow.processor.TaskPath) {
      // single input file allocate 2GB
      "2 GB"
    } else if (blastn_results instanceof nextflow.util.BlankSeparatedList) {
      // multiple input files
      // mem reqs = half of sum of GB file sizes plus 2GB wiggle room
      mem_reqs = Math.ceil(0.5 * (blastn_results.collect { it.size() }.sum()) / (1024**3)) + 2
      "${mem_reqs} GB" 
    } else {
      // not TaskPath or BlankSeparatedList, then default 2GB for memory
      "2 GB"
    }
  }
  conda 'conda-forge::python=3.10 conda-forge::biopython=1.80 conda-forge::openpyxl=3.1.0 conda-forge::pandas=1.5.3 conda-forge::rich=12.6.0 conda-forge::typer=0.7.0 conda-forge::xlsxwriter=3.0.8 conda-forge::polars=0.17.9 conda-forge::pyarrow=11.0.0'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
    }

    input:
    path excel_file
    path fasta_dir

    output:
    path "H1N1_samples.fasta", emit: h1
    path "H3N2_samples.fasta", emit: h3
    path "H5N1_samples.fasta", emit: h5
    path "VIC_samples.fasta", emit: vic
    path "YAM_samples.fasta", emit: yam

    script:
    """
    python3 sort_h_typing.py ${excel_file} ${fasta_dir}
    """
}


  
