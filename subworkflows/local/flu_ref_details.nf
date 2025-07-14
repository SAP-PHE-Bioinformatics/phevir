
include { BLAST_BLASTN as BLAST_BLASTN_IRMA                   } from '../../modules/local/blastn'
include { BLAST_BLASTN as BLAST_BLASTN_CONSENSUS              } from '../../modules/local/blastn'
include { BLAST_BLASTN as BLAST_BLASTN_CONSENSUS_REF_DB       } from '../../modules/local/blastn'




workflow REF_MAPPED_ANALYSIS{

take:
    blastn_irma
    meta
    ref_fasta

main:
    ch_versions = Channel.empty()

PULL_TOP_REF_ID(blastn_irma, meta.map{meta,file-> file}.first())
  ch_versions = ch_versions.mix(PULL_TOP_REF_ID.out.versions)

  PULL_TOP_REF_ID.out.accession_id
    .map { it[1] }
    .splitCsv(header: false, sep:",")
    // 0: sample_name, 1: segment, 2: ref_ncbi_accession_id, 3: ref_sequence_name
    .map{ [it[0], it[1], it[2]] }
    .combine(ch_filtered.map { [it[0].id, it[1]] }, by: 0)
    .set { ch_sample_segment } // ch_sample_segment: [sample_name, segment, id, reads]



  // Pull segment reference sequence for each sample
   SEQTK_SEQ(ch_sample_segment, ref_fasta.map{meta, file -> file}.first())
   ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions)
  
    // Map reads against segment reference sequences
    MINIMAP2(SEQTK_SEQ.out.sample_info)
    ch_versions = ch_versions.mix(MINIMAP2.out.versions)


    MOSDEPTH_GENOME(MINIMAP2.out.alignment)
    ch_versions = ch_versions.mix(MOSDEPTH_GENOME.out.versions)

    // Variants calling
        if (params.clair3_user_variant_model) {
        CLAIR3(
            MINIMAP2.out.alignment,
            ch_user_clair3_model
        )
        } else {
        CLAIR3(MINIMAP2.out.alignment, [])
        }
    ch_versions = ch_versions.mix(CLAIR3.out.versions)

    BCF_FILTER_CLAIR3(CLAIR3.out.vcf, params.major_allele_fraction)
    ch_versions = ch_versions.mix(BCF_FILTER_CLAIR3.out.versions)
    ch_vcf_filter = BCF_FILTER_CLAIR3.out.vcf


    VCF_FILTER_FRAMESHIFT(ch_vcf_filter)
    ch_versions = ch_versions.mix(VCF_FILTER_FRAMESHIFT.out.versions)

    BCFTOOLS_STATS(VCF_FILTER_FRAMESHIFT.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    VCF_FILTER_FRAMESHIFT.out.vcf
        .combine(MOSDEPTH_GENOME.out.bedgz, by: [0, 1, 2]) // combine channels based on sample_name, segment and accession_id
        .set { ch_bcf_consensus } // ch_bcf_consensus: [sample_name, segment, id, fasta, filt_vcf, mosdepth_per_base]

  
    // Generate consensus sequences
    BCF_CONSENSUS(ch_bcf_consensus, params.low_coverage)
    ch_versions = ch_versions.mix(BCF_CONSENSUS.out.versions)

    COVERAGE_PLOT(ch_bcf_consensus, params.low_coverage, ch_coverage_tsv)
    ch_versions = ch_versions.mix(COVERAGE_PLOT.out.versions)

    BCF_CONSENSUS.out.fasta
        .groupTuple(by: 0)
        .set { ch_final_consensus }

    CAT_CONSENSUS(ch_final_consensus)
    ch_versions = ch_versions.mix(CAT_CONSENSUS.out.versions)

    CAT_CONSENSUS.out.fasta
        .map { [[id:it[0]], it[1]] }
        .set { ch_cat_consensus }

    BLAST_BLASTN_CONSENSUS(ch_cat_consensus, BLAST_MAKEBLASTDB.out.db.first())
    ch_versions = ch_versions.mix(BLAST_BLASTN_CONSENSUS.out.versions)

    ch_blastn_consensus = BLAST_BLASTN_CONSENSUS.out.txt.collect({ it[1] })
    SUBTYPING_REPORT_BCF_CONSENSUS(
        PIGZ_META.out.file.map{meta,file-> file},
        ch_blastn_consensus,
        Channel.fromPath(params.input, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(SUBTYPING_REPORT_BCF_CONSENSUS.out.versions)
    subtyping_report = SUBTYPING_REPORT_BCF_CONSENSUS.out.report

    if (params.ref_db){
        BLAST_MAKEBLASTDB_REFDB(CHECK_REF_FASTA.out.fasta)
        ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB_REFDB.out.versions)

        BLAST_BLASTN_CONSENSUS_REF_DB(ch_cat_consensus, BLAST_MAKEBLASTDB_REFDB.out.db)
        ch_versions = ch_versions.mix(BLAST_BLASTN_CONSENSUS_REF_DB.out.versions)

        BLASTN_REPORT(BLAST_BLASTN_CONSENSUS_REF_DB.out.txt)
        ch_versions = ch_versions.mix(BLASTN_REPORT.out.versions)
    }

    datasets =Channel.of(
            [dataset: 'flu_h3n2_ha', reference: '', tag: ''],
            [dataset: 'flu_h1n1pdm_ha', reference: '', tag: ''],
            [dataset: 'community/moncla-lab/iav-h5/ha/all-clades', reference: '', tag: ''],
            [dataset: 'flu_vic_ha', reference: '', tag: ''],
            [dataset: 'flu_yam_ha', reference: '', tag: '']

        )
        .map { params -> tuple(params.dataset, params.reference, params.tag) }
    //nextclade download datasets
    NEXTCLADE_DATASETGET_FLU(datasets.map { it[0] }, datasets.map { it[1] }, datasets.map { it[2] })
    
    //Use channel ch_final_consensus to get just the HA and NA genes 
    //and then run nextclade on the tuple of the dataset and the consensus fasta for the HA and NA genes. 
    //ch_final_consensus.view()

    ch_ha_genes = ch_final_consensus
        .map { sample ->
            def sampleid = sample[0]
            def haFile = sample[1].find { it.getName().contains('HA') }
            [sampleid, haFile]
        }
        .filter { it[1] != null } // This ensures that only entries with a HA file are included  ha_genes.view()
    

    // Need to run each sample in ha_genes against each of the datasets in datasets
    ha_genes_dataset = ch_ha_genes.combine(NEXTCLADE_DATASETGET_FLU.out.dataset)
    //ha_genes_dataset.view()



    emit:

    ch_versions
    ch_coverage_tsv
    ch_final_consensus
    ch_cat_consensus
    subtyping_report
   }