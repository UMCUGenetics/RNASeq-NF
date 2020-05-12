include RSeQC from '../NextflowModules/RSeQC/3.0.1/RSeQC.nf' params( singleEnd:params.singleEnd)
include LcExtrap from '../NextflowModules/Preseq/2.0.3/LcExtrap.nf' params( optional:params.preseq.toolOptions)

workflow post_mapping_QC {
    take:
      bams_in
      bed_file
    main:
      RSeQC(bams_in, bed_file)
      LcExtrap(bams_in)
    emit:
      RSeQC.out, emit: rseqc_out
      LcExtrap.out, emit: preseq_out
}
