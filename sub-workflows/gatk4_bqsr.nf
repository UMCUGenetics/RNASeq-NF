include BaseRecalibrationTable from '../NextflowModules/GATK/4.1.3.0/BaseRecalibrationTable.nf' params(params)
include GatherBaseRecalibrationTables from '../NextflowModules/GATK/4.1.3.0/GatherBaseRecalibrationTables.nf' params(params)
include BaseRecalibration from '../NextflowModules/GATK/4.1.3.0/BaseRecalibration.nf' params(params)
include MergeBams from '../NextflowModules/Sambamba/0.6.8/MergeBams.nf' params(params)

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
