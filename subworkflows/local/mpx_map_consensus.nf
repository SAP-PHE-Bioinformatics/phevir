/*
============================================================================================================
    Assembly, Typing and Clade Variables Subworkflow Modules
============================================================================================================
*/
include { MINIMAP2_INDEX                      } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                      } from '../../modules/nf-core/minimap2/align/main'
include { IVAR_TRIM                           } from '../../modules/nf-core/ivar/trim/main'
include { SAMTOOLS_FAIDX                      } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_SORT                       } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                      } from '../../modules/nf-core/samtools/index/main'
include { IVAR_CONSENSUS                      } from '../../modules/nf-core/ivar/consensus/main'
include { BEDTOOLS_MASKFASTA                  } from '../../modules/nf-core/bedtools/maskfasta/main'
include { PANDEPTH                            } from '../../modules/local/pandepth.nf'

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

workflow MAP_CONSENSUS {
    take:
    clean_reads 
    ch_trim_ref     
    ch_reference           
    ch_primer                             
    ch_insert      

    main:
    ch_versions = Channel.empty()                        
                         


    //index reference
    MINIMAP2_INDEX(ch_trim_ref)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    //align reads to reference
    MINIMAP2_ALIGN(
        clean_reads,
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

    //rename bedfile chromosome name
    BEDTOOLS_MASKFASTA(
        params.bed_file,
        IVAR_CONSENSUS.out.fasta
    )

    PANDEPTH(
        ch_sorted_trimmed_bam.map{meta,bam,bai -> tuple([meta,bam])}
    )
    
    PANDEPTH.out.stats
        .groupTuple(by: [0][0])
        .map { id, stats_files ->
            def output_file = "${id}.stats.gz"
            """
            cat ${stats_files.collect { it[1] }.join(' ')} > ${output_file}
            """
            return tuple(id, file(output_file)) 
            }
        .set { grouped_stats }

    grouped_stats
    .collectFile(name: 'all_samples.stats.gz', storeDir: params.outdir) { id, stats_file ->
        stats_file
    }
    .set { all_samples_stats }



    emit:
    versions = ch_versions
    all_samples_stats
    sorted_trimmed_bam= ch_sorted_trimmed_bam
    fasta = BEDTOOLS_MASKFASTA.out.fasta
    pandepth_report = grouped_stats
    ref_index = SAMTOOLS_INDEX.out.bai

}