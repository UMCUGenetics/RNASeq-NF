include { MultiQC } from params.nextflowmodules_path+'/MultiQC/1.14/MultiQC.nf' params(optional: "--interactive --config $baseDir/assets/mqc_config.yaml")
include { RNASeqNFQC } from '../utils/RNASeqNFQC.nf' params(params)

workflow qc_report {
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
      //MultiQC report
      MultiQC(title, qc_files)
      //CustomQC report
      if (params.customQC) {
         RNASeqNFQC(title, Channel.fromPath(params.rmd_template), qc_files)
      } 
}

