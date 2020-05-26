include Quant from '../NextflowModules/Salmon/1.2.1/Quant.nf' params( singleEnd: params.singleEnd,
                                                                      stranded: params.stranded,
                                                                      unstranded: params.unstranded,
                                                                      revstranded: params.revstranded,
                                                                      saveUnaligned: params.saveUnaligned,
								      optional: params.options.Salmon_quant )
include QuantMerge from '../NextflowModules/Salmon/1.2.1/QuantMerge.nf' params( optional: params.options.Salmon_quantmerge  )
include MergeFastqLanes from '../NextflowModules/Utils/MergeFastqLanes.nf' params( params )
include Index as SalmonIndex from '../NextflowModules/Salmon/1.2.1/Index.nf' params( gencode: params.gencode,
                                                                                    optional: params.options.Salmon_index )
≈
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
      Quant ( MergeFastqLanes (fastq_files ), salmon_index.collect() )
      QuantMerge ( run_name, Quant.out.map { it[1] }.collect() )
      
    emit:
      quants_merged = QuantMerge.out
      logs =  Quant.out.quant_table 
}
