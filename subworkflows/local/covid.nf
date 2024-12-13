//modules
include { PANGOLIN } from '../../modules/nf-core/pangolin/main'
include { PANDEPTH } from '../../modules/local/pandepth.nf'
include { QUAST } from '../../modules/nf-core/quast/main'
include { ARTIC_MINION } from '../../modules/nf-core/artic/minion/main'
include { PANGO_COLLAPSE } from '../../modules/local/pangocollapse'

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
    //take the name of the primer bed file before the .bed as the metadata
    ch_primer_bed = Channel.of(['covid_primer', file(params.covid_primer_bed)])
    ch_reference = Channel.of(['covid_reference', file(params.covid_reference)])
    ch_covid_gff = Channel.of(['covid_gff', file(params.covid_gff)])
    ch_medaka_model = Channel.of(['medaka_model', file(params.medaka_model)]).view()

    ARTIC_MINION(
        READ_PREPROCESS.out.clean_reads,
        [],
        [],
        ch_reference.map { it[1] }.first(),
        ch_primer_bed.map { it[1] }.first(),
        ch_medaka_model.map{it[1]}.first(),
        '',
        params.scheme_name,
        params.scheme_version
    )
    ch_versions = ch_versions.mix(ARTIC_MINION.out.versions)

    PANDEPTH(
        ARTIC_MINION.out.bam
    )
    ch_versions = ch_versions.mix(PANDEPTH.out.versions) 

    
    PANDEPTH.out.stats
    .view() // Debugging step to inspect the channel with meta and file
    .map { meta, file -> file}
    .collectFile(name: 'COVID_pandepth_summary.tsv', storeDir: params.outdir, keepHeader: true)
    .set { pandepthSummaryFile }

pandepthSummaryFile.view() // Inspect the final result

    QUAST(
        ARTIC_MINION.out.fasta,
        ch_reference.first(),
        ch_covid_gff
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)

    QUAST.out.trans_tsv
    .map { it[1] }
    .collectFile(name: "${runID}_quast.tsv", storeDir: params.outdir, sort: true)
    .set { ch_quast_report }

    ARTIC_MINION.out.fasta
    .collectFile(name: 'covid_samples.fasta', storeDir: params.outdir) // Concatenate into one file
    .set { multiFasta }

    NEXTCLADE_SUB(
        channel.of([params.covid_nextclade_dataset_name]),
        multiFasta.map { [runID, it] }
    )
    ch_versions = ch_versions.mix(NEXTCLADE_SUB.out.versions)

    nextclade_summary = NEXTCLADE_SUB.out.report_tsv.collectFile(name: "${runID}_COVID_nextclade_summary.tsv", storeDir: params.outdir, sort: true)
    
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
    ch_pangolin_report = PANGOLIN.out.report.collectFile(name: "${runID}_pangolin.csv", storeDir: "${params.outdir}/pangolin", sort: true)
   
    ANNOTATION(
        ARTIC_MINION.out.fasta,
        params.covid_vadr_args,
        params.covid_vadr_trim_args
    )   
    ch_versions = ch_versions.mix(ANNOTATION.out.versions)


    // SUMMARY(
    // pandepthSummaryFile,
    // ch_quast_report,
    // PANGO_COLLAPSE.out.collapse_report,
    // nextclade_summary,
    // ch_pangolin_report     
    // )









    emit:
    ch_quast_report

}