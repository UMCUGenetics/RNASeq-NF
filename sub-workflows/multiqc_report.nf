include MultiQC from '../NextflowModules/MultiQC/1.8/MultiQC.nf' params(optional:params.multiqc.optional)

workflow multiqc_report {
    get:
      fastp_logs
      star_logs
      post_mapping_qc_logs
      htseq_logs 
      brecal_tables
    main: 
      qc_files = Channel.empty().mix( fastp_logs, star_logs, post_mapping_qc_logs, htseq_logs, brecal_tables).collect()
      qc_files.view()
      MultiQC(qc_files)
}
