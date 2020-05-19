include TrimGalore from '../NextflowModules/TrimGalore/0.6.5/TrimGalore.nf' params( optional: params.trimgalore.toolOptions, 
                                                                                   singleEnd: params.singleEnd )
include SortMeRNA from '../NextflowModules/SortMeRNA/4.2.0/SortMeRNA.nf' params( singleEnd:params.singleEnd )

workflow pre_processing {
    take:
      fastq_files
      sortmerna_fasta
    main:
      // Determine final fastqs files
      if ( params.runTrimGalore && params.runSortMeRna ) {
          TrimGalore(fastq_files) 
          SortMeRNA(TrimGalore.out.fastqs_trimmed, sortmerna_fasta.collect())
          final_fastqs = SortMeRNA.out.non_rRNA_fastqs
          trim_logs = TrimGalore.out.trimming_report
          fastqc_logs = TrimGalore.out.fastqc_report 
          srna_log = SortMeRNA.out.qc_report




      } else if ( params.runTrimGalore && !params.runSortMeRna ) {
          TrimGalore(fastq_files)
          final_fastqs = TrimGalore.out.fastqs_trimmed 
          trim_logs = TrimGalore.out.trimming_report
          fastqc_logs = TrimGalore.out.fastqc_report 
          srna_log = Channel.empty()


      } else if ( !params.runTrimGalore &&  params.runSortMeRna ) {
          SortMeRNA(fastq_files, sortmerna_fasta.collect() )
          final_fastqs = SortMeRNA.out.non_rRNA_fastqs 
          trim_logs = Channel.empty()
          fastqc_logs = Channel.empty()

      } else {
          final_fastqs = fastq_files
      }
      emit:
        processed_fastqs = final_fastqs
        srna_log = srna_log
        trim_logs = trim_logs
        fastqc_logs = fastqc_logs 
}
