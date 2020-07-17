include MultiQC from '../NextflowModules/MultiQC/1.9/MultiQC.nf' params( optional:params.options.MultiQC )
include RNASeqNFQC from '../utils/RNASeqNFQC.nf' params(params)
workflow multiqc_report {
    take:
      title
      fastqc_logs
      trim_logs
      star_logs
      post_mapping_qc_logs
      fc_logs
      salmon_logs
      sortmerna_logs
      flagstat_logs
     
      
    main: 
      qc_files = Channel.empty().mix( fastqc_logs, trim_logs, star_logs, post_mapping_qc_logs, flagstat_logs, fc_logs, salmon_logs, sortmerna_logs ).collect()
      MultiQC(title, qc_files)
      RNASeqNFQC(title, qc_files)
}
