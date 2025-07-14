process ZCAT_STATS {
    tag "$id"
    label 'process_low'
    container 'biocontainers/pigz:2.8'

    input:
    tuple val(id), path(files) 

    output:
    tuple val(id), path("*.stats"), emit: sample_stats

    
    
    script:
    """

    pigz -dc ${files[0]} | head -1 | sed 's/^#//' > ${id}.stats

   
    # Process all files, skipping headers for all except the first
    for file in ${files.join(' ')}; do
        pigz -dc "\$file" | awk -v sample="${id}" '
        BEGIN { OFS="\\t" }
        \$1 !~ /^#/ { print sample, \$0 } # Skip lines starting with "#" (header/footer)
        ' >> ${id}.stats
    done
    """
    
}

process MERGE_STATS {
    tag "${meta[0]}_${meta[1]}"
    label 'process_low'
    container  'quay.io/biocontainers/pandas:1.5.2'

    input:
    tuple val(meta), path(pandepth_file), path(coverage_file)

    output:
    path("${meta[0]}_${meta[1]}_merged.stats"), emit : merged_stats

    script:
    def prefix = "${meta[0]}_${meta[1]}"
    """
    merge_stats.py $pandepth_file $coverage_file ${meta[0]} ${meta[1]}    #pandepth_file, coverage_file, prefix, sample_name
    """
    }


process CONCAT_FILTER_CONSENSUS {
    tag { "${batch_id}-concat-filter" }
    label 'process_medium'

    input:
    tuple val(meta), path(summary_tsv) // Input summary TSV file

    output:
    path("${batch_id}_all.cat.consensus.fa") // Concatenated all fasta
    path("${batch_id}_passed_seqs.txt") // Passed sequences
    path("${batch_id}_at.cat.consensus.fa") // Austrakka concatenated fasta
    path("${batch_id}.gisaid.cat.consensus.fa") // GISAID concatenated fasta

    script:
    """
    # Concatenate all fasta files
    awk '/^>/ {gsub(/.consensus.fa(sta)?\$/, "", FILENAME); printf(">%s\\n", FILENAME); next;} {print}' 2_fasta/SA*.consensus.fasta > ${batch_id}_all.cat.consensus.fa
    sed -i 's,2_fasta/,,g' ${batch_id}_all.cat.consensus.fa

    # Filter sequences that pass criteria
    awk -F '\\t' '{if(\$4 >= 90 && \$2 >= 100 && \$1 !~/Sample/ && \$1 !~/NEG/ && \$1 !~/POS/) print \$1}' ${summary_tsv} > ${batch_id}_passed_seqs.txt

    # Create Austrakka concatenated fasta
    awk '{ if ((NR>1)&&(\$0~/^>/)) { printf("\\n%s", \$0); } else if (NR==1) { printf("%s", \$0); } else { printf("\\t%s", \$0); } }' ${batch_id}_all.cat.consensus.fa | grep -Ff ${batch_id}_passed_seqs.txt - | tr "\\t" "\\n" > ${batch_id}_at.cat.consensus.fa
    sed -i 's,2_fasta/,,g' ${batch_id}_at.cat.consensus.fa

    # Create GISAID concatenated fasta
    cat ${batch_id}_passed_seqs.txt | while read i; do
        awk '/^>/ {gsub(/.consensus.fa(sta)?\$/, "", FILENAME); printf(">%s\\n", "hCoV-19/Australia/" FILENAME "/2024"); next;} {print}' 2_fasta/"\$i".consensus.fasta >> /phe/viro/SARS-CoV-2/gisaid/${batch_id}.gisaid.cat.consensus.fa
    done
    sed -i 's,2_fasta/,,g' /phe/viro/SARS-CoV-2/gisaid/${batch_id}.gisaid.cat.consensus.fa
    """
}