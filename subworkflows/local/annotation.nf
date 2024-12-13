include { VADR; VADR_SUMMARIZE_ISSUES} from '../../modules/local/vadr'
include { PRE_TABLE2ASN; TABLE2ASN; POST_TABLE2ASN } from '../../modules/local/table2asn'


workflow ANNOTATION {
    take:
    assembly
    vadr_args
    vadr_trim_args
    
    main:
    ch_versions = Channel.empty()
    
    VADR(assembly, vadr_args, vadr_trim_args)
    ch_versions = ch_versions.mix(VADR.out.versions)
    VADR.out.feature_table
    .combine(VADR.out.pass_fasta, by: 0)
    .set { ch_pre_table2asn }
    
    VADR_SUMMARIZE_ISSUES(VADR.out.vadr_outdir.map { [it[1]] }.collect())
    ch_versions = ch_versions.mix(VADR.out.versions) 

    PRE_TABLE2ASN(ch_pre_table2asn)
    ch_versions = ch_versions.mix(PRE_TABLE2ASN.out.versions)
    TABLE2ASN(PRE_TABLE2ASN.out.table2asn_input)
    ch_versions = ch_versions.mix(TABLE2ASN.out.versions)
    POST_TABLE2ASN(TABLE2ASN.out.genbank)
    ch_versions = ch_versions.mix(POST_TABLE2ASN.out.versions)


    emit:
    genbank = POST_TABLE2ASN.out.genbank
    versions = ch_versions

}