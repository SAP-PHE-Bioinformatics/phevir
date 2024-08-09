#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIQC_TSV_FROM_LIST as READ_COUNT_FAIL_TSV        } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC_TSV_FROM_LIST as READ_COUNT_PASS_TSV        } from '../modules/local/multiqc_tsv_from_list'
include { PULL_TOP_REF_ID                                     } from '../modules/local/pull_top_ref_id'
include { IRMA                                                } from '../modules/local/irma'
include { SUBTYPING_REPORT as SUBTYPING_REPORT_IRMA_CONSENSUS } from '../modules/local/subtyping_report'
include { SUBTYPING_REPORT as SUBTYPING_REPORT_BCF_CONSENSUS  } from '../modules/local/subtyping_report'
include { COVERAGE_PLOT                                       } from '../modules/local/coverage_plot'
include { BLASTN_REPORT                                       } from '../modules/local/blastn_report'
include { VCF_FILTER_FRAMESHIFT                               } from '../modules/local/vcf_filter_frameshift'
include { MINIMAP2                                            } from '../modules/local/minimap2'
include { BCF_FILTER as BCF_FILTER_CLAIR3                     } from '../modules/local/bcftools'
include { BCF_FILTER as BCF_FILTER_MEDAKA                     } from '../modules/local/bcftools'
include { BCF_CONSENSUS; BCFTOOLS_STATS                       } from '../modules/local/bcftools'
include { CLAIR3                                              } from '../modules/local/clair3'
include { MOSDEPTH_GENOME                                     } from '../modules/local/mosdepth'
include { CAT_NANOPORE_FASTQ                                  } from '../modules/local/misc'
include { ZSTD_DECOMPRESS as ZSTD_DECOMPRESS_FASTA; ZSTD_DECOMPRESS as ZSTD_DECOMPRESS_CSV } from '../modules/local/zstd_decompress'
include { CAT_DB                                              } from '../modules/local/misc'
include { CAT_CONSENSUS                                       } from '../modules/local/misc'
include { SEQTK_SEQ                                           } from '../modules/local/seqtk_seq'
include { CHECK_SAMPLE_SHEET                                  } from '../modules/local/check_sample_sheet'
include { CHECK_REF_FASTA                                     } from '../modules/local/check_ref_fasta'

include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_NCBI         } from '../modules/local/blast_makeblastdb'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_REFDB        } from '../modules/local/blast_makeblastdb'
include { BLAST_BLASTN as BLAST_BLASTN_IRMA                   } from '../modules/local/blastn'
include { BLAST_BLASTN as BLAST_BLASTN_CONSENSUS              } from '../modules/local/blastn'
include { BLAST_BLASTN as BLAST_BLASTN_CONSENSUS_REF_DB       } from '../modules/local/blastn'
include { CUSTOM_DUMPSOFTWAREVERSIONS  as SOFTWARE_VERSIONS   } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { SORT_H_TYPING                                       } from '../modules/local/sort_ha_typing'
include { NEXTCLADE_DATASETGET                                 } from '../modules/local/nextclade_datasetget'
include { NEXTCLADE_RUN                                       } from '../modules/local/nextclade_run'
include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_NEXTCLADE          } from '../modules/local/multiqc_tsv_from_list'
include { MULTIQC                                             } from '../modules/local/multiqc'
include { NANOPLOT                                           } from '../modules/nf-core/nanoplot'
// include { BCFTOOLS_CSQ                                            } from '../modules/local/bcftools'

include { PIGZ as PIGZ_FASTA }  from '../modules/local/pigz'
include { PIGZ as PIGZ_META } from '../modules/local/pigz'
include { CREATE_REPORT } from '../modules/local/report'
include { CAT_FASTA } from '../modules/local/misc'
include { RESISTANCE } from '../modules/local/gisaid_resistance'
include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/kraken2/main'

//get last part of path as this is runID
runID = params.outdir.tokenize('/')[-1]
println "RunID: ${runID}"
def pass_sample_reads = [:]
def fail_sample_reads = [:]
ch_coverage_tsv = file(params.outdir + '/coverage.tsv')
ch_influenza_db_fasta = file(params.ncbi_influenza_fasta)
ch_influenza_metadata = file(params.ncbi_influenza_metadata)
if (params.clair3_user_variant_model) {
  ch_user_clair3_model = file(params.clair3_user_variant_model, checkIfExists: true)
}
def irma_module = 'FLU-minion'
if (params.irma_module) {
    irma_module = params.irma_module
}
def json_schema = "$projectDir/nextflow_schema.json"
def summary_params = NfcoreSchema.params_summary_map(workflow, params, json_schema)

workflow INFLUENZA {

    take:
    // TODO nf-core: edit input (take) channels
    ch_reads // channel: [ val(meta), [ fastq ] ]

    main:

    ch_for_multiqc = Channel.empty()
    ch_for_summary = Channel.empty()
    ch_versions    = Channel.empty()

  ch_input = CHECK_SAMPLE_SHEET(Channel.fromPath(params.input, checkIfExists: true))

  ch_input
    .map { [meta, fastq_1] }
    //count number of reads for sample
    .map { meta, fastq_1 -> [meta, fastq_1, fastq_1.countFastq()] }
    .set { ch_input_sorted }

  ch_input_sorted
    .branch { meta.id, fqgz, count  ->
      pass: count >= params.min_sample_reads
        pass_sample_reads[meta.id] = count
        return [ "$meta.id\t$count" ]
      fail: count < params.min_sample_reads
        fail_sample_reads[meta.id] = count
        return [ "$meta.id\t$count" ]
    }
    .set { ch_pass_fail_read_count }

  // Report samples which have reads count < min_sample_reads
  READ_COUNT_FAIL_TSV(
    ch_pass_fail_read_count.fail.collect(),
    ['Sample', 'Read count'],
    'fail_read_count_samples'
  )
  // Report samples which have reads count >= min_sample_reads
  READ_COUNT_PASS_TSV(
    ch_pass_fail_read_count.pass.collect(),
    ['Sample', 'Read count'],
    'pass_read_count_samples'
  )

  FASTP(
    ch_input_sorted,
    params.tfs_primers
    false,
    true,
    false
  )
  ch_versions = ch_versions.mix(FASTP.out.versions)

  NANOPLOT(FASTP.out.reads)
  ch_versions = ch_versions.mix(NANOPLOT.out.versions)
  // Keep samples which have reads count  > min_sample_reads for downstream analysis
  // Re-arrange channels to have meta map of information for sample
  ch_input_sorted
    .filter { it[-1] >= params.min_sample_reads }
    .map { meta, fastq, count -> [ [id: meta.id], fastq ] }
    .set { ch_reads }

  PIGZ_FASTA(file(params.ncbi_influenza_fasta), "influenza_ncbi.fasta")
  PIGZ_META(file(params.ncbi_influenza_metadata), "influenza_ncbi_metadata.csv")
  ch_versions = ch_versions.mix(PIGZ_FASTA.out.versions)

  ch_input_ref_db = PIGZ_FASTA.out.file

  if (params.ref_db){
    ch_ref_fasta = file(params.ref_db, type: 'file')
    CHECK_REF_FASTA(ch_ref_fasta)
    ch_versions = ch_versions.mix(CHECK_REF_FASTA.out.versions)
    CAT_DB(PIGZ_FASTA.out.file, CHECK_REF_FASTA.out.fasta)
    ch_input_ref_db = CAT_DB.out.fasta
  }

  BLAST_MAKEBLASTDB_NCBI(ch_input_ref_db)
  ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB_NCBI.out.versions)


  KRAKEN2_KRAKEN2(
        FASTP.out.reads.map{meta,fastqz -> [meta,fastqz]},
        params.kraken2_db,
        false,
        false
    )
  ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)
  ch_kraken_multiqc = KRAKEN2_KRAKEN2.out.report

  // IRMA to generate amended consensus sequences
  IRMA(FASTP.out.reads, irma_module)
  ch_versions = ch_versions.mix(IRMA.out.versions)
  // Find the top map sequences against ncbi database
  BLAST_BLASTN_IRMA(IRMA.out.majority_consensus, BLAST_MAKEBLASTDB_NCBI.out.db)
  ch_versions = ch_versions.mix(BLAST_BLASTN_IRMA.out.versions)

  ch_blast_irma = BLAST_BLASTN_IRMA.out.txt.collect({ it[1] })
  SUBTYPING_REPORT_IRMA_CONSENSUS(
    PIGZ_META.out.file,
    ch_blast_irma,
    CHECK_SAMPLE_SHEET.out
    )
  ch_versions = ch_versions.mix(SUBTYPING_REPORT_IRMA_CONSENSUS.out.versions)

  PULL_TOP_REF_ID(BLAST_BLASTN_IRMA.out.txt, PIGZ_META.out.file)
  ch_versions = ch_versions.mix(PULL_TOP_REF_ID.out.versions)

  PULL_TOP_REF_ID.out.accession_id
    .map { it[1] }
    .splitCsv(header: false, sep:",")
    // 0: sample_name, 1: segment, 2: ref_ncbi_accession_id, 3: ref_sequence_name
    .map { [id: it[0], segment: it[1], ref_id: it[2]] }
    .join(FASTP.out.reads, by: [0].id)
    // ch_sample_segment: [[sample_name, segment, ref_id], reads]
    .set { ch_sample_segment } 


  // Pull segment reference sequence for each sample
  SEQTK_SEQ(
            ch_sample_segment,
            ch_input_ref_db
)
  ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions)

  SEQTK_SEQ.out.sample_info
    .map{}
  
    MINIMAP2_ALIGN(
        FASTP.out.reads,
        ch_trim_ref.join(MINIMAP2_INDEX.out.index).first(),
        true,
        'bai',
        false,
        true
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    SAMTOOLS_FLAGSTAT(MINIMAP2_ALIGN.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)


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

    COVERAGE_PLOT(ch_bcf_consensus, params.low_coverage, ch_coverage_tsv)
    ch_versions = ch_versions.mix(COVERAGE_PLOT.out.versions)

    // Generate consensus sequences
    BCF_CONSENSUS(ch_bcf_consensus, params.low_coverage)
    ch_versions = ch_versions.mix(BCF_CONSENSUS.out.versions)

    BCF_CONSENSUS.out.fasta
        .groupTuple(by: 0)
        .set { ch_final_consensus }

    CAT_CONSENSUS(ch_final_consensus)
    ch_versions = ch_versions.mix(CAT_CONSENSUS.out.versions)

    CAT_CONSENSUS.out.fasta
        .map { [[id:it[0]], it[1]] }
        .set { ch_cat_consensus }

    BLAST_BLASTN_CONSENSUS(ch_cat_consensus, BLAST_MAKEBLASTDB_NCBI.out.db)
    ch_versions = ch_versions.mix(BLAST_BLASTN_CONSENSUS.out.versions)

    ch_blastn_consensus = BLAST_BLASTN_CONSENSUS.out.txt.collect({ it[1] })
    SUBTYPING_REPORT_BCF_CONSENSUS(
        PIGZ_META.out.file, 
        ch_blastn_consensus,
        CHECK_SAMPLE_SHEET.out
    )
    ch_versions = ch_versions.mix(SUBTYPING_REPORT_BCF_CONSENSUS.out.versions)

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
    NEXTCLADE_DATASETGET(datasets.map { it[0] }, datasets.map { it[1] }, datasets.map { it[2] })
    
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
    ha_genes_dataset = ch_ha_genes.combine(NEXTCLADE_DATASETGET.out.dataset)
    //ha_genes_dataset.view()
    
    VADR(ch_final_consensus)
    ch_versions = ch_versions.mix(VADR.out.versions)


    NEXTCLADE_RUN(ha_genes_dataset)
    
    NEXTCLADE_RUN.out.csv
                .map { 
                    meta, csv ->
                        def clade = Workflow.getNextcladeFieldMapFromCsv(csv)['clade']
                        return [ "$meta\t$clade" ]
                }
                .collect()                
                .map { 
                    tsv_data ->
                        def header = ['Sample', 'clade']
                        Workflow.multiqcTsvFromList(tsv_data, header)
                }
                .collectFile(name: "nextclade_run.tsv", storeDir: params.outdir)
                .set { ch_nextclade_multiqc }

    //cat irma consensus sequences

    //CAT_FASTA(runID, CAT_CONSENSUS.out.consensus_fasta)
    CAT_CONSENSUS.out.consensus_fasta
            .collectFile(name: "${runID}.fasta", storeDir: params.outdir)
            .set { consensus_for_run }

    //consensus_for_run.view()
    //uses the consensus fasta file to get the resistance mutations from FLUserver (GISAID)
    RESISTANCE(consensus_for_run)

    workflow_summary    = Schema.params_summary_multiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yaml")
    ch_kraken_multiqc.collect{it[1]}
    SOFTWARE_VERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

    MULTIQC(
        ch_multiqc_config,
        MINIMAP2.out.stats.collect().ifEmpty([]),
        MOSDEPTH_GENOME.out.mqc.collect().ifEmpty([]),
        BCFTOOLS_STATS.out.stats.collect().ifEmpty([]),
        ch_kraken_multiqc.collect{it[1]}.ifEmpty([]),
        SOFTWARE_VERSIONS.out.mqc_yml.collect(),
        ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
    )

    CREATE_REPORT(
        RESISTANCE.out.mutation_report,
        RESISTANCE.out.drug_sensitivity_report,
        SUBTYPING_REPORT_BCF_CONSENSUS.out.report,
        ch_nextclade_multiqc, 
        runID, 
        READ_COUNT_PASS_TSV.out.mqc_tsv,
        ch_coverage_tsv
    )
    
    }
