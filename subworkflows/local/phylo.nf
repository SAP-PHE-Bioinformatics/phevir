include { mafft }       from '../../modules/local/mafft'       addParams(params)
include { heatcluster } from '../../modules/local/heatcluster' addParams(params)
include { iqtree2 }     from '../../modules/local/iqtree'     addParams(params)
include { phytreeviz }  from '../../modules/local/phytreeviz'  addParams(params)
include { snpdists }    from '../../modules/local/snp-dists'   addParams(params)

workflow PHYLO {
  take:
    ch_fasta
    ch_reference_genome
    ch_prealigned
    
  main:
    ch_msa      = Channel.empty()
    ch_nwk      = Channel.empty()
    ch_versions = Channel.empty()
    ch_fasta.view()
    if ( params.msa == 'nextclade' ) {
      ch_msa = ch_prealigned.map { it -> tuple(it[0])}

    } 
    else if ( params.msa == 'mafft' ) {
      mafft(ch_fasta.map{it -> tuple(it[1])}.collect(), ch_reference_genome.map{ it -> tuple(it[1])})
      ch_msa = mafft.out.msa
      ch_versions = ch_versions.mix(mafft.out.versions)
    }
    
    iqtree2(ch_msa)
    phytreeviz(iqtree2.out.newick)
    snpdists(ch_msa)
    heatcluster(snpdists.out.matrix)
    ch_versions = ch_versions.mix(iqtree2.out.versions).mix(phytreeviz.out.versions).mix(snpdists.out.versions).mix(heatcluster.out.versions)

  emit:
    tree        = iqtree2.out.newick
    matrix      = snpdists.out.matrix
    msa         = ch_msa
    for_multiqc = phytreeviz.out.for_multiqc.mix(heatcluster.out.for_multiqc)
    versions    = ch_versions
}