/*
====================================================================================================
    Nextclade Dataset and Analysis Subworkflow Modules
====================================================================================================
*/

include { NEXTCLADE_DATASETGET          } from '../../modules/local/nextclade/nextclade_datasetget.nf'
include { NEXTCLADE                     } from '../../modules/local/nextclade_run.nf'
include { NEXTCLADE_CAT                  } from '../../modules/local/nextclade_cat.nf'
/*
====================================================================================================
    Run Nextclade Dataset and Analysis Subworkflow
====================================================================================================
*/

workflow NEXTCLADE_SUB {
    take:
    dataset
    run_samples

    main:

    //debug input
    //run_samples.view()

    ch_versions              = Channel.empty()
    ch_nextclade_report      = Channel.empty()
    ch_prealigned         = Channel.empty()
    ch_nextclade_run_input   = Channel.empty()

    if (params.skip_nextclade) return // conditional check on param.skip_nextclade. If true, subworkflow will not execute.

    
    NEXTCLADE_DATASETGET(dataset)
    ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)

    input= run_samples.join(NEXTCLADE_DATASETGET.out.dataset_fetch)

    NEXTCLADE(input.map { it[1] }, input.map { [it[0], it[2]] })
    ch_prealigned.mix(NEXTCLADE.out.prealigned)
    ch_nextclade_report = NEXTCLADE.out.nextclade_file

    // NEXTCLADE_CAT(NEXTCLADE.out.tsv)
    // sorted_tsv_files = NEXTCLADE_CAT.out.nextclade_cat_tsv


    emit:
    fasta_aligned          = NEXTCLADE.out.prealigned
    report_tsv           = NEXTCLADE.out.tsv
    report_csv           = NEXTCLADE.out.nextclade_file
    nextclade_report       = ch_nextclade_report
    versions               = ch_versions
}