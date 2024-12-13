#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PULL_TOP_REF_ID                                    } from '../../modules/local/pull_top_ref_id'
include { IRMA                                                } from '../../modules/local/irma'
include { SUBTYPING_REPORT as SUBTYPING_REPORT_IRMA_CONSENSUS } from '../../modules/local/subtyping_report'
include { SUBTYPING_REPORT as SUBTYPING_REPORT_BCF_CONSENSUS  } from '../../modules/local/subtyping_report'
include { BLASTN_REPORT                                       } from '../../modules/local/blastn_report'
include { VCF_FILTER_FRAMESHIFT                               } from '../../modules/local/vcf_filter_frameshift'
include { MINIMAP2                                       } from '../../modules/local/minimap2_custom'
include { BCF_FILTER as BCF_FILTER_CLAIR3                     } from '../../modules/local/bcftools'
include { BCF_CONSENSUS; BCFTOOLS_STATS                       } from '../../modules/local/bcftools'
include { CLAIR3                                              } from '../../modules/local/clair3'
include { MOSDEPTH_GENOME                                     } from '../../modules/local/mosdepth'
include { CAT_DB                                              } from '../../modules/local/misc'
include { CAT_CONSENSUS                                       } from '../../modules/local/misc'
include { SEQTK_SEQ                                           } from '../../modules/local/seqtk_seq'
include { CHECK_REF_FASTA                                     } from '../../modules/local/check_ref_fasta'

include { BLAST_MAKEBLASTDB      } from '../../modules/local/blast_makeblastdb'
include { BLAST_BLASTN as BLAST_BLASTN_IRMA                   } from '../../modules/local/blastn'
include { SORT_H_TYPING                                       } from '../../modules/local/sort_ha_typing'
//include { MULTIQC_TSV_FROM_LIST as MULTIQC_TSV_NEXTCLADE          } from '../../modules/local/multiqc_tsv_from_list'
include { MULTIQC                                             } from '../../modules/local/multiqc'
include { FASTP                                               } from '../../modules/nf-core/fastp'
include { PIGZ_UNCOMPRESS as PIGZ_FASTA; PIGZ_UNCOMPRESS as PIGZ_META } from '../../modules/nf-core/pigz/uncompress'
include { SUMMARY } from '../../modules/local/summary_flu'
include { CAT_FASTA } from '../../modules/local/misc'
include { RESISTANCE } from '../../modules/local/gisaid_resistance'
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/main'

include { KRAKENTOPMATCH } from '../../modules/local/kraken2parser'

include { ASSEMBLY_QC } from '../../subworkflows/local/flu_assembly_qc'
include { READ_PREPROCESS } from '../../subworkflows/local/read_filt'
include { NEXTCLADE_SUB } from '../../subworkflows/local/nextclade_sub'
include { ANNOTATION } from '../../subworkflows/local/annotation'
include { REF_MAPPED_ANALYSIS } from '../../subworkflows/local/flu_ref_details'
include { NEXTCLADE_DATASETGET } from '../../modules/local/nextclade/nextclade_datasetget'
include { NEXTCLADE as NEXTCLADE_RUN } from '../../modules/local/nextclade_run'

//get last part of path as this is runID


workflow INFLUENZA {

    take:

    ch_reads // channel: [ val(meta), [ fastq ] ]

    main:

    runID = params.run_id
    println "RunID: ${runID}"
    ch_coverage_tsv = file(params.outdir + '/coverage.tsv')
    ch_influenza_db_fasta = Channel.of(['NCBI_db', file(params.ncbi_influenza_fasta)]).view()
    ch_influenza_metadata = Channel.of(['NCBI_meta', file(params.ncbi_influenza_metadata)])
    if (params.clair3_user_variant_model) {
      ch_user_clair3_model = file(params.clair3_user_variant_model, checkIfExists: true)
    }
    
    

    ch_for_multiqc = Channel.empty()
    ch_for_summary = Channel.empty()
    ch_versions    = Channel.empty()
    ch_filtered    = Channel.empty()


  PIGZ_FASTA(ch_influenza_db_fasta)
  PIGZ_META(ch_influenza_metadata)
  ch_versions = ch_versions.mix(PIGZ_FASTA.out.versions)

  ch_input_ref_db = PIGZ_FASTA.out.file

  BLAST_MAKEBLASTDB(ch_input_ref_db)
  ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

  READ_PREPROCESS(ch_reads, params.kraken2_db)
  ch_versions = ch_versions.mix(READ_PREPROCESS.out.versions)


  ASSEMBLY_QC(READ_PREPROCESS.out.clean_reads)
  ch_versions = ch_versions.mix(ASSEMBLY_QC.out.versions)


  // Find the top map sequences against ncbi database
  BLAST_BLASTN_IRMA(ASSEMBLY_QC.out.assembly, BLAST_MAKEBLASTDB.out.db.first())
  ch_versions = ch_versions.mix(BLAST_BLASTN_IRMA.out.versions)

  ch_blast_irma = BLAST_BLASTN_IRMA.out.txt.collect({ it[1] })
  SUBTYPING_REPORT_IRMA_CONSENSUS(
    PIGZ_META.out.file.map{meta,file-> file},
    ch_blast_irma,
    Channel.fromPath(params.input, checkIfExists: true)
    )
  ch_versions = ch_versions.mix(SUBTYPING_REPORT_IRMA_CONSENSUS.out.versions)


  if(params.map_to_reference){
   REF_MAPPED_ANALYSIS(ASSEMBLY_QC.out.assembly, PIGZ_META.out.file, PIGZ_FASTA.out.file)
   ch_versions = ch_versions.mix(REF_MAPPED_ANALYSIS.out.versions)
  }
  
    ch_datasets =Channel.of(
            ['flu_h3n2_ha', 'flu_h3n2_ha'],
            ['flu_h1n1pdm_ha', 'flu_h1n1pdm_ha'],
            ['flu_h5_all', 'community/moncla-lab/iav-h5/ha/all-clades'],
            ['flu_vic_ha', 'flu_vic_ha'],
            ['flu_yam_ha', 'flu_yam_ha']

        )

ASSEMBLY_QC.out.HA
    .map { sample ->
        def meta = sample[0]
        def file = sample[1]
        return [meta.dataset, file] // Pair dataset with its FASTA file
    }
    .groupTuple() // Group FASTA files by dataset
    .set { groupedFastaByDataset }


  groupedFastaByDataset
    .map { dataset, fastaFiles ->
        // Create a concatenated FASTA file
        def baseName = new File(dataset).name
        def concatFastaPath = "grouped_${baseName}.fasta"
        new File(concatFastaPath).withWriter { writer ->
            fastaFiles.each { file ->
                 file.eachLine { line ->
                    writer << line + "\n" // Append each line with a newline
                }
            }
        }
        return [dataset, file(concatFastaPath)] // Wrap the path with file()
    }
    .set { nextcladeInputs }


    NEXTCLADE_SUB(
        nextcladeInputs.map { it[0] }, // Dataset channel
        nextcladeInputs  // FASTA channel
    )

   ANNOTATION(ASSEMBLY_QC.out.assembly, params.flu_vadr_args, params.flu_vadr_trim_args)


    //CAT_FASTA(runID, CAT_CONSENSUS.out.consensus_fasta)
    ASSEMBLY_QC.out.assembly
            .map { it[1] }
            .collectFile(name: "${runID}.fasta", storeDir: params.outdir, sort: true)
            .set { consensus_for_run }


    //uses the consensus fasta file to get the resistance mutations from FLUserver (GISAID)
    RESISTANCE(consensus_for_run)

    nextclade_summary = NEXTCLADE_SUB.out.report_tsv.collectFile(name: "${runID}_nextclade_summary.tsv", storeDir: params.outdir, sort: true, keepHeader: true, skip: 1)

  ch_multiqc = Channel.empty()
  
    SUMMARY(
        ASSEMBLY_QC.out.stats,
        ASSEMBLY_QC.out.irma_consensus_qc,
        ASSEMBLY_QC.out.typing_report,
        RESISTANCE.out.mutation_report,
        RESISTANCE.out.drug_sensitivity_report,
        SUBTYPING_REPORT_IRMA_CONSENSUS.out.report,
        nextclade_summary,
        runID, 
        READ_PREPROCESS.out.kraken2_reportsheet_tsv
    )
    
    

    emit:
    multiqc = ch_multiqc
    //report = SUMMARY.out.report
    versions = ch_versions
    }
