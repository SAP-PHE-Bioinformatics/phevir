/*
========================================================================================================
    Preprocessing Read QC Subworkflow Modules
========================================================================================================
*/

include { NCBI_SRA_HUMAN_SCRUBBER              } from '../../modules/local/ncbi_sra_human_scrubber.nf'

include { KRAKEN2_KRAKEN2                      } from '../../modules/nf-core/kraken2/main.nf'
include { KRAKEN2REPORT_SUMMARY                } from '../../modules/local/kraken2report_summary.nf'
include { KRAKEN2_REPORTSHEET                  } from '../../modules/local/kraken2_reportsheet.nf'
//include { QC_REPORT                            } from '../../modules/local/qc_report.nf'
include { NANOSTAT                             } from '../../modules/nf-core/nanostat/main.nf'
include { CHOPPER                              } from '../../modules/nf-core/chopper/main.nf'


/*
========================================================================================================
    Run Preprocessing Read QC Subworkflow
========================================================================================================
*/

workflow READ_PREPROCESS {
    take:
    reads
    db // params.krakendb

    main:
    ch_versions                = Channel.empty()
    ch_kraken2reportsheet      = Channel.empty()
    ch_kraken2_reportsheet_tsv = Channel.empty()

    
    NCBI_SRA_HUMAN_SCRUBBER(reads)
    ch_versions = ch_versions.mix(NCBI_SRA_HUMAN_SCRUBBER.out.versions)


    CHOPPER(NCBI_SRA_HUMAN_SCRUBBER.out.reads)
    ch_filtered = CHOPPER.out.fastq

    NANOSTAT(ch_filtered)
    ch_versions = ch_versions.mix(NANOSTAT.out.versions)
    

      KRAKEN2_KRAKEN2(
        ch_filtered.map{meta,fastqz -> [meta,fastqz]},
        db,
        false,
        false
    )
   ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)
   
   ch_kraken2report_summary_input = KRAKEN2_KRAKEN2.out.report
   KRAKEN2REPORT_SUMMARY(ch_kraken2report_summary_input)
   ch_kraken2reportsheet = KRAKEN2REPORT_SUMMARY.out.kraken_lines.collect()
   KRAKEN2_REPORTSHEET(ch_kraken2reportsheet)
   // Populate ch_kraken2_reportsheet_tsv with actual data
   ch_kraken2_reportsheet_tsv = KRAKEN2_REPORTSHEET.out.kraken2_reportsheet_tsv


    emit:
    clean_reads                = CHOPPER.out.fastq
    stats                      = NANOSTAT.out.txt.collect()
    versions                   = ch_versions
    kraken2_reportsheet_tsv    = KRAKEN2_REPORTSHEET.out.kraken2_reportsheet_tsv
}