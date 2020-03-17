include MultiQC from '../NextflowModules/MultiQC/1.8/MultiQC.nf' params(params)

workflow multiqc_report {
    get:
      fastp_logs
      star_logs
      post_mapping_qc_logs
      htseq_logs
    main: 
      qc_files = Channel.empty().mix( fastp_logs, star_logs, post_mapping_qc_logs, htseq_logs).collect()
      qc_files.view()
      MultiQC(qc_files)
}
