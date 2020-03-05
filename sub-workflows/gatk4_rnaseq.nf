include SplitIntervals from '../NextflowModules/GATK/4.1.3.0/SplitIntervals.nf' params(params)
include SplitNCigarReads from '../NextflowModules/GATK/4.1.3.0/SplitNCigarReads.nf' params(params)
include BaseRecalibrationTable from '../NextflowModules/GATK/4.1.3.0/BaseRecalibrationTable.nf' params(params)
include GatherBaseRecalibrationTables from '../NextflowModules/GATK/4.1.3.0/GatherBaseRecalibrationTables.nf' params(params)
include BaseRecalibration from '../NextflowModules/GATK/4.1.3.0/BaseRecalibration.nf' params(params)
include MergeBams from '../NextflowModules/Sambamba/0.6.8/MergeBams.nf' params(params)
include HaplotypeCaller from '../NextflowModules/GATK/4.1.3.0/HaplotypeCaller.nf' params(params)
include CombineGVCFs from '../NextflowModules/GATK/4.1.3.0/CombineGVCFs.nf' params(params)

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
      vcf_combine_per_interval = HaplotypeCaller.out
         .groupTuple(by:[1])
         .map{
           sample_ids, interval, gvcfs, idxs, interval_files ->
           [run_id, interval, gvcfs, idxs, interval_files[0]]
         } 
       /* Combine GVCFs per interval */
      CombineGVCFs(
       vcf_combine_per_interval
      )

    emit:
      bams_recal = MergeBams.out
      vcf = CombineGVCFs.out
}
