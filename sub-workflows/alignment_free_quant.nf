include { Quant } from params.nextflowmodules_path+'/Salmon/1.9.0/Quant.nf' params( single_end: params.single_end,
                                                                      stranded: params.stranded,
                                                                      unstranded: params.unstranded,
                                                                      revstranded: params.revstranded,
                                                                      saveUnaligned: params.saveUnaligned,
								      optional: params.options.Salmon_quant )
include { QuantMerge } from params.nextflowmodules_path+'/Salmon/1.9.0/QuantMerge.nf' params( optional: params.options.Salmon_quantmerge  )
include { Index as SalmonIndex } from params.nextflowmodules_path+'/Salmon/1.9.0/Index.nf' params( gencode: params.gencode,
                                                                                    optional: params.options.Salmon_index )

workflow alignment_free_quant {
    take:
      fastq_files
      run_name

    main:
      if ( params.salmon_index ) {
          salmon_index = Channel
              .fromPath(params.salmon_index, checkIfExists: true)
              .ifEmpty { exit 1, "Transcripts fasta not found: ${params.salmon_index}"}
      } else if ( !params.salmon_index ) {
          transcripts_fasta = Channel
              .fromPath(params.transcripts_fasta, checkIfExists: true)
              .ifEmpty { exit 1, "Fasta file not found: ${params.transcripts_fasta}"}
          SalmonIndex ( transcripts_fasta )
          salmon_index = SalmonIndex.out.salmon_index
      }
      Quant ( fastq_files.map { [it[0],it[2]] }, salmon_index.collect() )
      QuantMerge ( run_name, Quant.out.map { it[1] }.collect() )
      
    emit:
      quants_merged = QuantMerge.out
      logs =  Quant.out.quant_table.map { it[1] } 
}
