/*
============================================================================================================
    Assembly, Typing and Clade Variables Subworkflow Modules
============================================================================================================
*/

include { IRMA                                 } from '../../modules/local/irma.nf'
include { IRMA_CONSENSUS_QC                    } from '../../modules/local/irma_consensus_qc.nf'
include { IRMA_REPORT       } from '../../modules/local/irma_report.nf'
include { IRMA_SEGMENT_COVERAGE                } from '../../modules/local/segment_coverage.nf'
include { MERGE_COVERAGE_RESULTS               } from '../../modules/local/merge_coverage_results.nf'
include { ABRICATE_FLU                         } from '../../modules/local/abricate_flu.nf'
include { ABRICATE_REPORT                 } from '../../modules/local/abricate_report.nf'
include { NEXTCLADE_SORT                  } from '../../modules/local/nextclade_sort.nf'
include { PANDEPTH                             } from '../../modules/local/pandepth.nf'
include { ZCAT_STATS                          } from '../../modules/local/misc_tasks.nf'
include { MERGE_STATS                          } from '../../modules/local/misc_tasks.nf'

/*
============================================================================================================
    Assembly, Typing and Clade Variables Subworkflow Params Setup
============================================================================================================
*/

def irma_module = ''
if (params.irma_module) {
    irma_module = params.irma_module
}

/*
============================================================================================================
    Run Assembly, Typing, and Clade Variables Subworkflow
============================================================================================================
*/

workflow ASSEMBLY_QC {
    take:
    clean_reads 

    main:
    ch_versions                        = Channel.empty()
    ch_assembly                        = Channel.empty()
    ch_HA                              = Channel.empty()
    ch_NA                              = Channel.empty()

    irma_module = Channel.of('FLU-minion')
    IRMA(clean_reads, irma_module.first())

    ch_assembly = IRMA.out.assembly

    ch_versions = ch_versions.mix(IRMA.out.versions)

    ch_HA = IRMA.out.HA
    ch_NA = IRMA.out.NA

    IRMA_CONSENSUS_QC(IRMA.out.assembly)
    ch_majority_consensus = IRMA.out.majority_consensus
    irma_consensus_qc_files = IRMA_CONSENSUS_QC.out.irma_consensus_qc

    ch_irma_consensus_qc_results = irma_consensus_qc_files
        .unique { meta, file_path -> meta.id }  // Use unique to remove duplicates, 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collect all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def qc_header = list[0].split("\n")[0]
            def qc_contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != qc_header }
            return ([qc_header] + qc_contentWithoutHeaders).join("\n")
        }

    IRMA_REPORT(ch_irma_consensus_qc_results)
    irma_consensus_qc_tsv = IRMA_REPORT.out.irma_consensus_qc_tsv.collectFile(name: "irma_consensus_qc.tsv", storeDir: params.outdir)

    IRMA.out.irma_fasta
        .flatMap { meta, fasta_file_paths ->
        // Collect each bam file and map it with the updated metadata
        fasta_file_paths.collect { fasta_file ->
            def file_name = fasta_file.name

            // Extract segment information from the file name
            def segment = file_name.replaceFirst(/A_/, '').replace('.fasta', '')

            // Add the segment to the meta information
            def updated_meta = meta + [segment: segment]

            // Return each [meta, fasta_file] as separate elements
            return [updated_meta, fasta_file]
        }
    }
   .set{ fasta_files_individual }

    IRMA.out.irma_bam
    .flatMap { meta, bam_file_paths ->
        // Collect each bam file and map it with the updated metadata
        bam_file_paths.collect { bam_file ->
            def file_name = bam_file.name

            // Extract segment information from the file name
            def segment = file_name.replaceFirst(/A_/, '').replace('.bam', '')

            // Add the segment to the meta information
            def updated_meta = meta + [segment: segment]

            // Return each [meta, bam_file] as separate elements
            return [updated_meta, bam_file]
        }
    }
    .set { bam_with_updated_meta }

    PANDEPTH(bam_with_updated_meta)
   

    IRMA_SEGMENT_COVERAGE(fasta_files_individual)
    irma_seg_cov_results_files = IRMA_SEGMENT_COVERAGE.out.cov_results

    ch_pandepth_mapped = PANDEPTH.out.stats
    .map { meta, file -> [[meta.id, meta.segment], file] }

    // Prepare the second channel similarly
    ch_coverage_mapped = IRMA_SEGMENT_COVERAGE.out.cov_results
    .map { meta, file -> [[meta.id, meta.segment], file] }
   
   // Join the channels by the keys (meta.id and meta.segment)
    ch_combined = ch_pandepth_mapped
    .join(ch_coverage_mapped, by: 0)


    MERGE_STATS(ch_combined)
    
    all_sample_stats=MERGE_STATS.out.merged_stats.collectFile(name: "all_samples_merged.stats", storeDir: params.outdir, keepHeader: true, sort: true, skip: 1)

    ABRICATE_FLU(IRMA.out.assembly)
    ch_versions = ch_versions.mix(ABRICATE_FLU.out.versions)

    ch_irma_abricate_report_input = IRMA.out.tsv.join(ABRICATE_FLU.out.tsv)

    ABRICATE_REPORT(ch_irma_abricate_report_input)
    tsv_files = ABRICATE_REPORT.out.tsv_combined

    ch_combined_results = tsv_files
        .unique { meta, file_path -> meta.id }  // Use unique to remove duplicates, 'id' is the unique key in meta
        .map { meta, file_path -> file_path.text }  // Convert each file to its textual content
        .flatten()  // Flatten the channel to process each line individually
        .filter { line -> line && line.trim() != '' }  // Filter out null or empty lines
        .collect()  // Collect all the lines into a list
        .map { list ->
            // Include the header only once at the start of the combined file
            def header = list[0].split("\n")[0]
            def contentWithoutHeaders = list*.split("\n").flatten().unique().findAll { it != header }
            return ([header] + contentWithoutHeaders).join("\n")
        }


    typing_report_tsv = ABRICATE_REPORT.out.tsv_combined.map { meta, file -> file }.collectFile(name: 'flu_typing_report.tsv', storeDir: params.outdir, sort: true, keepHeader: true, skip: 1)

    ch_nextclade_for_sort = ABRICATE_REPORT.out.tsv_combined

    NEXTCLADE_SORT(ch_nextclade_for_sort)

    // Join IRMA.out.HA and ASSEMBLY_QC.out.dataset channels on the 'id' field
    IRMA.out.HA.join(NEXTCLADE_SORT.out.dataset, by: [0][0] )
    .map { meta_irma, ha_file, dataset_file ->
        def subtype = dataset_file.text.trim() // Extract subtype from the dataset file
        def updated_meta = meta_irma + [dataset: subtype] // Append subtype to the meta
        return [updated_meta, ha_file] // Return updated meta and HA_file
    }
    .set { irma_meta_HA }



    emit:
    HA                              = irma_meta_HA
    NA                              = IRMA.out.NA
    typing_report              = typing_report_tsv
    stats                           = all_sample_stats
    irma_consensus_qc          = irma_consensus_qc_tsv
    assembly                        = ch_assembly
    versions                        = ch_versions

}