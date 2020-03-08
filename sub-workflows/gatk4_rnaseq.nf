include SplitIntervals from '../NextflowModules/GATK/4.1.3.0/SplitIntervals.nf' params(params)
include SplitNCigarReads from '../NextflowModules/GATK/4.1.3.0/SplitNCigarReads.nf' params(params)
include BaseRecalibrationTable from '../NextflowModules/GATK/4.1.3.0/BaseRecalibrationTable.nf' params(params)
include GatherBaseRecalibrationTables from '../NextflowModules/GATK/4.1.3.0/GatherBaseRecalibrationTables.nf' params(params)
include BaseRecalibration from '../NextflowModules/GATK/4.1.3.0/BaseRecalibration.nf' params(params)
include MergeBams from '../NextflowModules/Sambamba/0.6.8/MergeBams.nf' params(params)
include HaplotypeCaller from '../NextflowModules/GATK/4.1.3.0/HaplotypeCaller.nf' params(params)
include MergeVCFs from '../NextflowModules/GATK/4.1.3.0/MergeVCFs.nf' params(params)
include VariantFiltration from '../NextflowModules/GATK/4.1.3.0/VariantFiltration.nf' params(params)

workflow gatk4_rnaseq {
    get:
      bam_dedup
      genome_fasta
      genome_index
      genome_dict
    main:
      SplitNCigarReads(bam_dedup, genome_fasta.collect(), genome_index.collect(), genome_dict.collect())
      SplitIntervals( 'no-break', Channel.fromPath( params.scatter_interval_list))
      BaseRecalibrationTable(SplitNCigarReads.out.combine(SplitIntervals.out.flatten()))
      GatherBaseRecalibrationTables(BaseRecalibrationTable.out.groupTuple())
      //Perform BQSR
      BaseRecalibration(
        SplitNCigarReads.out
           .combine(GatherBaseRecalibrationTables.out, by:0)
           .combine(SplitIntervals.out.flatten())
      )
      //Merge recalibrated bams
      MergeBams(
        BaseRecalibration.out
           .groupTuple()
           .map{ [it[0],it[2],it[3]] }
      )
      HaplotypeCaller( MergeBams.out.combine(SplitIntervals.out.flatten()))
      MergeVCFs(
        HaplotypeCaller.out.groupTuple(by:[0]).map{
          sample_id, intervals, gvcfs, idxs, interval_files ->
          [sample_id, gvcfs, idxs]
        }
      )
      VariantFiltration( MergeVCFs.out.map{
          sample_id, vcfs, idxs -> [sample_id, "filtered", "RNA", vcfs, idxs] }
      )
    emit:
      bams_recal = MergeBams.out
      vcf = MergeVCFs.out
      filter = VariantFiltration.out 
}
