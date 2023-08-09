include { EdgerNormalize } from params.nextflowmodules_path+'/edgeR/3.40.0/normalize.nf' params( params )
include { Deseq2Normalize } from params.nextflowmodules_path+'/DESeq2/1.38.0/deseq2Normalize.nf' params( params )
include { FeatureCounts } from params.nextflowmodules_path+'/Subread/2.0.3/FeatureCounts.nf' params( optional: params.options.FeatureCounts,
                                                                                       biotypeQC:params.biotypeQC,
                                                                                       single_end: params.single_end,
                                                                                       stranded: params.stranded,
                                                                                       unstranded: params.unstranded,
                                                                                       revstranded: params.revstranded,
                                                                                       fc_group_features: params.fc_group_features,
                                                                                       fc_count_type: params.fc_count_type,
                                                                                       fc_group_features_type: params.fc_group_features_type,
                                                                                       fc_extra_attributes: params.fc_extra_attributes,
                                                                                       gencode: params.gencode)

workflow alignment_based_quant {
    take:
        run_name
        bam_files
        genome_gtf

    main:
        FeatureCounts( run_name, bam_files.collect(), genome_gtf.collect() )
        if ( params.normalize_counts ) {
            EdgerNormalize (run_name, FeatureCounts.out.count_table )
            Deseq2Normalize (run_name, FeatureCounts.out.count_table )
        }
      
    emit:
        fc_read_counts =  FeatureCounts.out.count_table
        fc_summary = FeatureCounts.out.count_summary
}
