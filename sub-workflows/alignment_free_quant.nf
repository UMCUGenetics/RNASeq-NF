include Count from '../NextflowModules/HTSeq/0.11.3/Count.nf' params(hts_count_type:params.hts_count_type,
                                                                    hts_group_features:params.hts_group_features,
                                                                    optional:params.count.toolOptions, 
                                                                    singleEnd:params.singleEnd, 
                                                                    stranded:params.stranded, 
                                                                    unstranded:params.unstranded, 
                                                                    revstranded:params.revstranded)

include Quant from '../NextflowModules/Salmon/1.2.1/Quant.nf' params(singleEnd: params.singleEnd,
                                                                     stranded: params.stranded,
                                                                     unstranded: params.unstranded,
                                                                     revstranded: params.revstranded,
                                                                     saveUnaligned: params.saveUnaligned,
								                                     optional: params.salmon_quant.toolOptions )
include QuantMerge from '../NextflowModules/Salmon/1.2.1/QuantMerge.nf' params( optional: params.salmon_quantmerge.toolOptions )
include MergeFastqLanes from '../NextflowModules/Utils/MergeFastqLanes.nf' params(params)


workflow alignment_free_quant {
    take:
      fastq_files
      salmon_index
      run_name
    main:
      Quant ( MergeFastqLanes (fastq_files ), salmon_index.collect() )
      QuantMerge ( run_name, Quant.out.map { it[1] }.collect() )
    emit:
      quants_merged = QuantMerge.out
      logs =  Quant.out.quant_table 
}