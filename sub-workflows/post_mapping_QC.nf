include RSeQC from '../NextflowModules/RSeQC/3.0.1/RSeQC.nf' params( singleEnd:params.singleEnd)
include LCExtrap from '../NextflowModules/Preseq/2.0.3/LCExtrap.nf' params( optional:params.preseq.toolOptions)

workflow post_mapping_QC {
    take:
      bams_in
      bed_file
    main:
      RSeQC(bams_in, bed_file)
      LCExtrap(bams_in)
    emit:
      RSeQC.out
      LCExtrap.out
}
