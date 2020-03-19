include HaplotypeCaller from '../NextflowModules/GATK/4.1.3.0/HaplotypeCaller.nf' params(params)
include VariantFiltration from '../NextflowModules/GATK/4.1.3.0/VariantFiltration.nf' params(params)
include MergeVCFs from '../NextflowModules/GATK/4.1.3.0/MergeVCFs.nf' params(params)

workflow gatk4_hc {
    get:
      bam
      scatter_intervals
    main:
      HaplotypeCaller(bam.combine(scatter_intervals))
      MergeVCFs(
        HaplotypeCaller.out.groupTuple(by:[0]).map{
          sample_id, intervals, gvcfs, idxs, interval_files ->
          [sample_id, gvcfs, idxs]
        }
      )
      VariantFiltration( MergeVCFs.out.map{
          sample_id, vcfs, idxs -> [sample_id, "all", "RNA", vcfs, idxs] }
      )
    emit:
      VariantFiltration.out 
}
