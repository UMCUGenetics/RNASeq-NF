include BaseRecalibrationTable from '../NextflowModules/GATK/4.1.3.0/BaseRecalibrationTable.nf' params(mem:params.baserecalibrator.mem,
												       optional:params.baserecalibrator.toolOptions,
												       genome_known_sites:params.genome_known_sites,
                                                                                                       genome_fasta:params.genome_fasta)
include GatherBaseRecalibrationTables from '../NextflowModules/GATK/4.1.3.0/GatherBaseRecalibrationTables.nf' params(mem:params.gatherbaserecalibrator.mem)
include BaseRecalibration from '../NextflowModules/GATK/4.1.3.0/BaseRecalibration.nf' params(mem:params.applybqsr.mem,
                                                                                             genome_fasta: params.genome_fasta)
include MergeBams from '../NextflowModules/Sambamba/0.6.8/MergeBams.nf' params(mem:params.mergebams.mem,
                                                                               optional:params.mergebams.toolOptions)

workflow gatk4_bqsr {
    get:
      bam
      scatter_intervals
    main:
      BaseRecalibrationTable(bam.combine(scatter_intervals))
      GatherBaseRecalibrationTables(BaseRecalibrationTable.out.groupTuple())
      //Perform BQSR
      BaseRecalibration(
        bam
          .combine(GatherBaseRecalibrationTables.out, by:0)
          .combine(scatter_intervals)
      )
      //Merge recalibrated bams
      MergeBams(
        BaseRecalibration.out
           .groupTuple()
           .map{ [it[0],it[2],it[3]] }
      )
    emit:
      MergeBams.out	
      GatherBaseRecalibrationTables.out
}
