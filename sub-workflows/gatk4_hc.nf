include HaplotypeCaller from '../NextflowModules/GATK/4.1.3.0/HaplotypeCaller.nf' params (mem:params.haplotypecaller.mem,
                                                                                          genome_fasta:params.genome_fasta,
                                                                                          optional:params.haplotypecaller.toolOptions )
include VariantFiltration from '../NextflowModules/GATK/4.1.3.0/VariantFiltration.nf' params(mem:params.variantfiltration.mem,,
                                                                                             genome_fasta:params.genome_fasta,
                                                                                             optional: params.variantfiltration.toolOptions )
include MergeVCFs as MergeVCF from '../NextflowModules/GATK/4.1.3.0/MergeVCFs.nf' params( mem: params.mergevcf.mem )


workflow gatk4_hc {
    take:
      bam
      scatter_intervals
      run_id
    main:
      HaplotypeCaller(bam.combine(scatter_intervals))
      MergeVCF(
        HaplotypeCaller.out.groupTuple(by:[0]).map{
          sample_id, intervals, gvcfs, idxs, interval_files ->
          [sample_id, gvcfs, idxs]
        }
      )
      //Filter raw vcf files
      VariantFiltration( MergeVCF.out.map{
        sample_id, vcfs, idxs -> [sample_id, run_id, "RNA", vcfs, idxs] }
      )
     
    emit:
      VariantFiltration.out
}
