//modules
include { PANGOLIN } from '../../modules/nf-core/pangolin/main'
include { PANDEPTH } from '../../modules/local/pandepth.nf'
include { QUAST } from '../../modules/nf-core/quast/main'
include { ARTIC_MINION } from '../../modules/nf-core/artic/minion/main'
include { PANGO_COLLAPSE } from '../../modules/local/pangocollapse'
include { SUMMARY_COV } from '../../modules/local/summary_covid'
include { SAMTOOLS_DEPTH_COV } from '../../modules/local/samtools_depth_cov'
include { DATAMASH } from '../../modules/local/datamash'


//subworkflows
include { NEXTCLADE_SUB } from './nextclade_sub'
include { READ_PREPROCESS } from './read_filt'
include { ANNOTATION } from './annotation'
// include ( SUMMARY ) from '../../modules/local/summary_covid'

nextflow.enable.dsl = 2

workflow COVID {

    take: 
    ch_reads

    main:

    ch_versions = Channel.empty()
    READ_PREPROCESS(ch_reads, params.kraken2_db)
    ch_versions = ch_versions.mix(READ_PREPROCESS.out.versions)

    runID= params.run_id
    amplicon = "COV"
    //take the name of the primer bed file before the .bed as the metadata
    ch_primer_bed = Channel.of(['covid_primer', file(params.covid_primer_bed)])
    ch_reference = Channel.of(['covid_reference', file(params.covid_reference)])
    ch_covid_gff = Channel.of(['covid_gff', file(params.covid_gff)])
    ch_medaka_model = Channel.of(['medaka_model', file(params.medaka_model)])
    ch_scheme = Channel.of(params.scheme_name)
    ch_scheme_version = Channel.of(params.scheme_version)

    ARTIC_MINION(
        READ_PREPROCESS.out.clean_reads,
        [],
        [],
        ch_reference.map { it[1] }.first(),
        ch_primer_bed.map { it[1] }.first(),
        ch_medaka_model.map { it[1] }.first(),
        ch_scheme.first(),
        ch_scheme_version.first()
    )
    ch_versions = ch_versions.mix(ARTIC_MINION.out.versions)

    PANDEPTH(
        ARTIC_MINION.out.bam
    )
    ch_versions = ch_versions.mix(PANDEPTH.out.versions) 

    
    // PANDEPTH.out.stats
    // //.view() // Debugging step to inspect the channel with meta and file
    // .map { meta, file -> file}
    // .collectFile(name: 'COVID_pandepth_summary.tsv', storeDir: "${params.outdir}/COVID", keepHeader: true)
    // .set { pandepthSummaryFile }
   PANDEPTH.out.stats
    .map { meta, file ->
        // Read the original file
        def lines = file.text.readLines()

        // Add Sample column to the header and prepend meta.id to valid data lines
        def updatedLines = lines.collect { line ->
            if (line.startsWith('#Chr')) {
                // Add "Sample" column if not already present in the header
                line.contains("Sample") ? line : "Sample\t${line}"
            } else if (line.startsWith('##')) {
                null // Drop comment lines starting with ##
            } else if (!line.startsWith(meta.id)) {
                // Prepend Sample ID only if not already added
                "${meta.id}\t${line}"
            } else {
                line // Keep the line as is
            }
        }.findAll { it != null } // Filter out null lines (dropped ## lines)

        // Write the modified content back to the same file
        file.text = updatedLines.join('\n') + '\n'

        // Return the modified file
        return file
    }
    .collectFile(name: 'COVID_pandepth_summary.tsv', storeDir: "${params.outdir}/${runID}_${amplicon}", keepHeader: true)
    .set { pandepthSummaryFile }

//pandepthSummaryFile.view() // Inspect the final result

    SAMTOOLS_DEPTH_COV (
        ARTIC_MINION.out.bam_primertrimmed
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH_COV.out.versions) 
    

    DATAMASH (
        SAMTOOLS_DEPTH_COV.out.tsv
    )
    ch_versions = ch_versions.mix(DATAMASH.out.versions) 
    
    QUAST(
        ARTIC_MINION.out.fasta,
        ch_reference.first(),
        ch_covid_gff.first()
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)

    QUAST.out.trans_tsv
    .map { it[1] }
    .collectFile(name: "${runID}_quast.tsv", storeDir: "${params.outdir}/${runID}_${amplicon}", sort: true)
    .set { ch_quast_report }

    ARTIC_MINION.out.fasta
    .map { it[1] }
    .collectFile(name: "${runID}_covid.fasta", storeDir: "${params.outdir}/${runID}_${amplicon}") // Concatenate into one file
    .set { multiFasta }

    NEXTCLADE_SUB(
        channel.of([params.covid_nextclade_dataset_name]),
        multiFasta.map { [params.covid_nextclade_dataset_name, it] }
    )
    ch_versions = ch_versions.mix(NEXTCLADE_SUB.out.versions)

    nextclade_summary = NEXTCLADE_SUB.out.report_tsv.collectFile(name: "${runID}_COVID_nextclade_summary.tsv", storeDir: "${params.outdir}/${runID}_${amplicon}", sort: true)
    
    PANGO_COLLAPSE(
        NEXTCLADE_SUB.out.report_tsv,
        params.collapse_file,
        runID
    )
    ch_versions = ch_versions.mix(PANGO_COLLAPSE.out.versions)

    PANGOLIN(
        multiFasta.map { [runID, it] }
    )
    ch_versions = ch_versions.mix(PANGOLIN.out.versions)
    ch_pangolin_report = PANGOLIN.out.report.collectFile(name: "${runID}_pangolin.csv", storeDir: "${params.outdir}/${runID}_COV/pangolin", sort: true)
   
    ANNOTATION(
        ARTIC_MINION.out.fasta,
        params.covid_vadr_args,
        params.covid_vadr_trim_args
    )   
    ch_versions = ch_versions.mix(ANNOTATION.out.versions)


    SUMMARY_COV(
    pandepthSummaryFile,
    ch_pangolin_report,   
    nextclade_summary,
    PANGO_COLLAPSE.out.collapsed,  
    ARTIC_MINION.out.fasta.map { it[1] }.collect(),
    runID
    )


    emit:
    ch_quast_report

}