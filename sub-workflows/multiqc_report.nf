include MultiQC from '../NextflowModules/MultiQC/1.8/MultiQC.nf' params(optional:params.multiqc.optional)

workflow multiqc_report {
    take:
      fastqc_logs
      trim_logs
      sortmerna_logs
      star_logs
      post_mapping_qc_logs
      htseq_logs 
      fc_logs
      salmon_logs
    main: 
      qc_files = Channel.empty().mix( fastqc_logs, trim_logs, sortmerna_logs, star_logs, post_mapping_qc_logs, htseq_logs, fc_logs, salmon_logs ).collect()
      MultiQC( qc_files )
}
