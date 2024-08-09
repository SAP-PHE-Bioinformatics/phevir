// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join


include { SAMTOOLS_SORT      } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX }  from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'
include { NANOPLOT            } from '../../modules/nf-core/nanoplot/main'
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
include { VADR               } from '../../modules/local/vadr'
include { NEXTCLADE_DATASETGET } from '../../modules/local/nextclade_datasetget'
include { NEXTCLADE          } from '../../modules/local/nextclade_run'
include { SUMMARY } from '../../modules/local/summary'
include { PHYLO } from './phylo.nf'
include { CHOPPER } from '../../modules/nf-core/chopper/main'


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

    //QC and filtering
    NANOPLOT(ch_reads)
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)
    
    CHOPPER(ch_reads)
    ch_versions = ch_versions.mix(CHOPPER.out.versions)

    //index reference
    MINIMAP2_INDEX(ch_trim_ref)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    //align reads to reference
    MINIMAP2_ALIGN(
        CHOPPER.out.fastq,
        ch_trim_ref.join(MINIMAP2_INDEX.out.index).first(),
        true,
        'bai',
        false,
        true
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    
    //trim reads
    IVAR_TRIM(
        MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index, by: [0][0]),
        ch_primer.map{it -> it[1]}.first()
    )
    ch_versions = ch_versions.mix(IVAR_TRIM.out.versions)

    //index masked reference
    SAMTOOLS_FAIDX ( ch_reference )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    //sort trimmed reads
    SAMTOOLS_SORT(
        IVAR_TRIM.out.bam,
        ch_reference.join(SAMTOOLS_FAIDX.out.fai).first()
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_sorted_trimmed_bam = SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0][0])
            .map (meta, bam, bai) -> {
             [ meta, bam, bai  ]
        }
    //ch_sorted_trimmed_bam.view()
    IVAR_CONSENSUS(
        ch_sorted_trimmed_bam.map{meta,bam,bai -> tuple([meta,bam])}.filter{meta, bam -> meta.species != 'NEG'},
        ch_reference.map{it -> it[1]}.first(),
        false
    )
    ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions)

    //qc

    KRAKEN2_KRAKEN2(
        ch_reads,
        params.kraken2_db,
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

    ACI(
        ch_sorted_trimmed_bam,
        ch_insert.map{it -> it[1]}.first()
    )
    ch_versions = ch_versions.mix(ACI.out.versions)

    SAMTOOLS_STATS(
        ch_sorted_trimmed_bam,
        ch_reference.first()
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    // SAMTOOLS_COVERAGE(
    //     ch_sorted_trimmed_bam,
    //     ch_reference.first(),
    //     SAMTOOLS_FAIDX.out.fai.first()   
    // )
    // ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)
    QUAST(
        IVAR_CONSENSUS.out.fasta,
        ch_reference.first(),
        tuple('MPOX', params.gff)
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)


    SAMTOOLS_DEPTH(
        ch_sorted_trimmed_bam.map{meta,bam,bai -> tuple([meta,bam])},
        ch_insert.first()
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)
    IVAR_VARIANTS(
        ch_sorted_trimmed_bam.map{meta, bam, bai -> tuple([meta,bam])}.filter{meta, bam -> meta.species != 'NEG'},
        ch_reference.map{it -> it[1]}.first(),
        SAMTOOLS_FAIDX.out.fai.map{it -> it[1]}.first(),       
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
        ch_sorted_trimmed_bam.combine(ch_primer),
    )   
    ch_versions = ch_versions.mix(SAMTOOLS_AMPLICONSTATS.out.versions)

    PLOT(
        SAMTOOLS_AMPLICONSTATS.out.samtools_ampliconstats_files,
    )
    ch_versions = ch_versions.mix(PLOT.out.versions.first())

    FASTQC(
        ch_reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    SAMTOOLS_FLAGSTAT(
        ch_sorted_trimmed_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    // SAMTOOLS_COVERAGE.out.coverage
    //   .collectFile(name: "samtools_coverage_summary.tsv",
    //     keepHeader: true,
    //     storeDir: "${params.outdir}/samtools_coverage")
    //   .set { samtools_coverage_file }

    ACI.out.cov
      .collectFile(name: "aci_coverage_summary.csv",
        keepHeader: true,
        storeDir: "${params.outdir}/aci")
      .set { aci_coverage_file }

    VADR(
        IVAR_CONSENSUS.out.fasta.map{it -> it[1]}.collect()
    )
    ch_versions = ch_versions.mix(VADR.out.versions)
    NEXTCLADE_DATASETGET()
    ch_versions = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
    NEXTCLADE(
        ch_for_tree
        .mix(IVAR_CONSENSUS.out.fasta)
        .map{meta,fasta -> fasta}
        .collect()
        , NEXTCLADE_DATASETGET.out.dataset
    )
    ch_versions = ch_versions.mix(NEXTCLADE.out.versions)
    ch_for_multiqc = ch_for_multiqc.mix(NEXTCLADE.out.nextclade_file)
    ch_for_summary.map{ it -> it[1]}.collect()
    ch_for_summary = ch_for_summary
        .mix(aci_coverage_file.map{it -> it[1]})
        .mix(SAMTOOLS_STATS.out.stats.map{it -> it[1]})
        //.mix(SAMTOOLS_COVERAGE.out.coverage.map{it -> it[1]})
        .mix(SAMTOOLS_DEPTH.out.tsv.map{it -> it[1]})
        .mix(IVAR_VARIANTS.out.tsv.map{it -> it[1]})
        .mix(SAMTOOLS_AMPLICONSTATS.out.samtools_ampliconstats_files.map{it -> it[1]})
        .mix(VADR.out.vadr_file)
        .mix(NEXTCLADE.out.nextclade_file)
        .mix(QUAST.out.tsv.map{it -> it[1]})

    //ch_for_summary.view()
    ch_for_multiqc = ch_for_multiqc.mix(SAMTOOLS_FLAGSTAT.out.flagstat.map{it -> it[1]}).mix(ACI.out.for_multiqc)

    

    PHYLO(
        ch_for_tree.concat(IVAR_CONSENSUS.out.fasta),
        ch_reference,
        NEXTCLADE.out.prealigned
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

