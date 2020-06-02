include MultiQC from '../NextflowModules/MultiQC/1.8/MultiQC.nf' params( optional:params.options.MultiQC )

workflow multiqc_report {
    take:
      title
      fastqc_logs
      trim_logs
      sortmerna_logs
      star_logs
      post_mapping_qc_logs
      fc_logs
     
      
    main: 
      qc_files = Channel.empty().mix( fastqc_logs, trim_logs, sortmerna_logs, star_logs, post_mapping_qc_logs, fc_logs ).collect()
      MultiQC(title, qc_files)
}
