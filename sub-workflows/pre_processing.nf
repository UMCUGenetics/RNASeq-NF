include TrimGalore from '../NextflowModules/TrimGalore/0.6.5/TrimGalore.nf' params( optional: params.options.TrimGalore, 
                                                                                   single_end: params.single_end )
include SortMeRNA from '../NextflowModules/SortMeRNA/4.3.3/SortMeRNA.nf' params( single_end: params.single_end )
include FastQC from '../NextflowModules/FastQC/0.11.8/FastQC.nf' params( optional: params.options.FastQC) 

workflow pre_processing {
    take:
      fastq_files
    main:
      //Empty log channels
      srna_logs = Channel.empty()
      trim_logs = Channel.empty()
      fastqc_logs = Channel.empty()
      //
      if (params.runFastQC) {
         FastQC(fastq_files)
         fastqc_logs = FastQC.out.report
      }
      if ( params.runSortMeRNA) {
           rRNA_database = file(params.rRNA_database_manifest)
           //if (rRNA_database.isEmpty()) {exit 1, "File ${rRNA_database.getName()} is empty!"}
           sortmerna_fasta = Channel
                .from( rRNA_database.readLines() )
                .map { row -> file(row) }
      }
      // Determine final fastqs files
      if ( params.runTrimGalore && params.runSortMeRNA ) {
          TrimGalore(fastq_files) 
          SortMeRNA(TrimGalore.out.fastqs_trimmed, sortmerna_fasta.collect())
          final_fastqs = SortMeRNA.out.non_rRNA_fastqs
          trim_logs = TrimGalore.out.trimming_report
          fastqc_logs = TrimGalore.out.fastqc_report 
          srna_logs = SortMeRNA.out.qc_report

      } else if ( params.runTrimGalore && !params.runSortMeRNA ) {
          TrimGalore(fastq_files)
          final_fastqs = TrimGalore.out.fastqs_trimmed 
          trim_logs = TrimGalore.out.trimming_report
          fastqc_logs = TrimGalore.out.fastqc_report 
         
      } else if ( !params.runTrimGalore &&  params.runSortMeRNA ) {
          SortMeRNA(fastq_files, sortmerna_fasta.collect() )
          final_fastqs = SortMeRNA.out.non_rRNA_fastqs 
          srna_logs = SortMeRNA.out.qc_report

      } else {
          final_fastqs = fastq_files
      }
      emit:
        processed_fastqs = final_fastqs
        srna_logs = srna_logs
        trim_logs = trim_logs
        fastqc_logs = fastqc_logs 
}
