include { RNASeqNFDGE } from '../utils/RNASeq_DGE.nf' params( params )    

workflow generate_dge_report {
    take:
        run_name
        count_file
        md_file
        comparisons_file
        dge_template

    main:
        RNASeqNFDGE( run_name, count_file, md_file, comparisons_file, dge_template )
      
    emit:
        dge_report =  RNASeqNFDGE.out.dge_report
}
