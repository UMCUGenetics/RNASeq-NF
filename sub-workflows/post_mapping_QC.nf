include RSeQC from '../../NextflowModules/RSeQC/2.6.1/RSeQC.nf' params(params)
include Lc_extrap from '../../NextflowModules/Preseq/2.0.3/Lc_extrap.nf' params(params)

workflow post_mapping_QC {
    get:
      bams_in
      bed_file
    main:
      RSeQC(bams_in, bed_file)
      Lc_extrap(bams_in)
    emit:
      RSeQC.out
      Lc_extrap.out
}
