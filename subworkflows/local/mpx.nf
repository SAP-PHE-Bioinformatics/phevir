// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join


include { SAMTOOLS_SORT      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX }  from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'
include { MINIMAP2_ALIGN      } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX      } from '../../modules/nf-core/minimap2/index/main'
include { IVAR_TRIM           } from '../../modules/nf-core/ivar/trim/main'
include { IVAR_CONSENSUS      } from '../../modules/nf-core/ivar/consensus/main'
include { IVAR_VARIANTS       } from '../../modules/nf-core/ivar/variants/main'
include { SAMTOOLS_STATS      } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_COVERAGE   } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_DEPTH      } from '../../modules/nf-core/samtools/depth/main'
include { QUAST }     from '../../modules/nf-core/quast'
include { SAMTOOLS_AMPLICONSTATS ; PLOT } from '../../modules/local/samtools_ampliconstat'
include { ACI                } from '../../modules/local/aci'
include { KRAKEN2_KRAKEN2           } from '../../modules/nf-core/kraken2/main'
include { FASTQC             } from '../../modules/nf-core/fastqc/main'
include { SUMMARY } from '../../modules/local/summary'
include { PHYLO } from './phylo.nf'
include { BEDTOOLS_MASKFASTA } from '../../modules/nf-core/bedtools/maskfasta/main'

include { READ_PREPROCESS } from '../../subworkflows/local/read_filt'
include { MAP_CONSENSUS } from '../../subworkflows/local/mpx_map_consensus'
include { PANDEPTH } from '../../modules/local/pandepth'
include { ANNOTATION } from './annotation'
include { NEXTCLADE_SUB } from './nextclade_sub'

workflow MPX {

    take:
    // TODO nf-core: edit input (take) channels
    ch_reads // channel: [ val(meta), [ fastq ] ]

    main:

    ch_for_multiqc = Channel.empty()
    ch_for_summary = Channel.empty()
    ch_versions    = Channel.empty()


    //steup and admin tasks 

    ch_reference = channel.of(tuple('MPOX', params.mpx_reference))
    ch_primer = channel.of(tuple('MPOX_primer', params.primer))
    ch_insert = channel.of(tuple('MPOX_primer', params.insert))
    ch_for_tree = channel.fromPath("${params.mpx_consensus_sequences}/*", type: 'file', checkIfExists: true)
        .map{ it -> tuple(tuple(id:"${it.baseName}_illumina"), it) }
    ch_trim_ref = Channel.fromPath("${params.trim_reference}/*", type: 'file', checkIfExists: true)
        .map{ it -> tuple(tuple(id:"${it.baseName}"), it) }

    READ_PREPROCESS(ch_reads, params.kraken2_db)
    ch_versions = ch_versions.mix(READ_PREPROCESS.out.versions)

    MAP_CONSENSUS(ch_reads, ch_trim_ref, ch_reference, ch_primer, ch_insert)
    ch_versions = ch_versions.mix(MAP_CONSENSUS.out.versions)

    ACI(
        MAP_CONSENSUS.out.sorted_trimmed_bam,
        ch_insert.map{it -> it[1]}.first()
    )
    ch_versions = ch_versions.mix(ACI.out.versions)

    SAMTOOLS_STATS(
        MAP_CONSENSUS.out.sorted_trimmed_bam,
        ch_reference.first()
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    QUAST(
        MAP_CONSENSUS.out.fasta,
        ch_trim_ref.first(),
        tuple('MPOX', params.gff)
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)

    IVAR_VARIANTS(
        MAP_CONSENSUS.out.sorted_trimmed_bam.map{meta, bam, bai -> tuple([meta,bam])}.filter{meta, bam -> meta.species != 'NEG'},
        ch_reference.map{it -> it[1]}.first(),
        MAP_CONSENSUS.out.ref_index.map{it -> it[1]}.first(),       
        params.gff,
        'true'
    )
    ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions)

    // IVAR_VARIANTS.out.vcf
    //     .join(ch_sorted_trimmed_bam)
    //     .combine(ch_reference)
    //     .set { for_igv_reports }
    // IGV_REPORTS(for_igv_reports)
    // ch_versions= ch_versions.mix(IGV_REPORTS.out.versions)

    SAMTOOLS_AMPLICONSTATS(
        MAP_CONSENSUS.out.sorted_trimmed_bam.combine(ch_primer),
    )   
    ch_versions = ch_versions.mix(SAMTOOLS_AMPLICONSTATS.out.versions)

    PLOT(
        SAMTOOLS_AMPLICONSTATS.out.samtools_ampliconstats_files,
    )
    ch_versions = ch_versions.mix(PLOT.out.versions.first())

    SAMTOOLS_FLAGSTAT(
        MAP_CONSENSUS.out.sorted_trimmed_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    ACI.out.cov
      .collectFile(name: "aci_coverage_summary.csv",
        keepHeader: true,
        storeDir: "${params.outdir}/aci")
      .set { aci_coverage_file }

    // VADR(MAP_CONSENSUS.out.fasta)
    // ch_versions = ch_versions.mix(VADR.out.versions)

    // VADR.out.feature_table
    // .combine(VADR_FLU.out.pass_fasta, by: 0)
    // .set { ch_pre_table2asn }
    
    // VADR_SUMMARIZE_ISSUES(VADR.out.vadr_outdir.map { [it[1]] }.collect())
    // ch_versions = ch_versions.mix(VADR.out.versions)
    ANNOTATION(
        MAP_CONSENSUS.out.fasta,
        params.mpx_vadr_args,
        params.mpx_vadr_trim_args
    ) 


    NEXTCLADE_SUB(
        channel.of([params.mpx_nextclade_dataset]),
        //ch_for_tree.mix(
        MAP_CONSENSUS.out.fasta
    )

    ch_for_summary.map{ it -> it[1]}.collect()
    // ch_for_summary = ch_for_summary
    //     .mix(aci_coverage_file.map{it -> it[1]})
    //     .mix(SAMTOOLS_STATS.out.stats.map{it -> it[1]})
    //     .mix(MAP_CONSENSUS.out.pandepth_report.map{it -> it[1]})
    //     .mix(IVAR_VARIANTS.out.tsv.map{it -> it[1]})
    //     .mix(SAMTOOLS_AMPLICONSTATS.out.samtools_ampliconstats_files.map{it -> it[1]})
    //     // .mix(VADR.out.vadr_file)
    //     .mix(NEXTCLADE.out.nextclade_file)
    //     .mix(QUAST.out.tsv.map{it -> it[1]})

    //ch_for_summary.view()
    ch_for_multiqc = ch_for_multiqc.mix(SAMTOOLS_FLAGSTAT.out.flagstat.map{it -> it[1]}).mix(ACI.out.for_multiqc)

    PHYLO(
        ch_for_tree.concat(MAP_CONSENSUS.out.fasta),
        ch_reference,
        NEXTCLADE_SUB.out.fasta_aligned
    )
    ch_versions = ch_versions.mix(PHYLO.out.versions)
    ch_for_multiqc = ch_for_multiqc.mix(PHYLO.out.for_multiqc)

    tree = PHYLO.out.tree
    alignment = PHYLO.out.msa
    matrix = PHYLO.out.matrix
    //     .mix(NEXTCLADE.out.for_multiqc)
    //     .set(multiqc)
    
    emit:
    // bam      = ch_sorted_trimmed_bam      
    // consensus = IVAR_CONSENSUS.out.fasta         
    summary  = ch_for_summary
    multiqc = ch_for_multiqc
    // tree     = tree
    // alignment = alignment
    // matrix   = matrix

    

    
    versions = ch_versions                     // channel: [ versions.yml ]
}

