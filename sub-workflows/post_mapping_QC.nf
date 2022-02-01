include RSeQC from '../NextflowModules/RSeQC/3.0.1/RSeQC.nf' params( single_end:params.single_end)
include RSeQC_TIN from '../NextflowModules/RSeQC/3.0.1/RSeQC.nf' params( params )
include LCExtrap from '../NextflowModules/Preseq/2.0.3/LCExtrap.nf' params( optional:params.options.Preseq )
include GtfToGenePred from '../NextflowModules/UCSC/377/GtfToGenePred.nf' params( params )
include GenePredToBed from '../NextflowModules/UCSC/377/GenePredToBed.nf' params( params )

include SkipLowReadBam from '../utils/SkipLowReadBam.nf' params( params )

workflow post_mapping_QC {
    take:
      bams_in
      genome_gtf

    main:
        rseqc_logs = Channel.empty()
        rseqc_tin_logs = Channel.empty()
        preseq_lce_logs = Channel.empty()

      //def bam_subset = bams_in.map { sample_id, rg_id, bam, bai -> [sample_id, bam, bai] }
      if (params.genome_bed ) {
          //Create bed12 index file
          genome_bed = Channel
              .fromPath(params.genome_bed, checkIfExists: true)
              .ifEmpty { exit 1, "Genome bed12 file not found: ${params.genome_bed}"}
      } else if ( !params.genome_bed ) {
          GtfToGenePred ( genome_gtf)
          GenePredToBed ( GtfToGenePred.out.genome_genepred )
          genome_bed = GenePredToBed.out.genome_bed12
      }
        //while RSeQC does complete without issues it causes problems with multiQC crashing in case a bam file without reads is supplied
      SkipLowReadBam(bams_in)

      //RSeQC(bams_in, genome_bed.collect())
      RSeQC(SkipLowReadBam.out.bam.map { sample_id, rg_id, bam, bai -> [sample_id, bam, bai] }, genome_bed.collect())
      
      if (params.runRSeQC_TIN) {
          //RSeQC_TIN(bams_in, genome_bed.collect()) 
          RSeQC_TIN(SkipLowReadBam.out.bam.map { sample_id, rg_id, bam, bai -> [sample_id, bam, bai] }, genome_bed.collect()) 
          rseqc_tin_logs = RSeQC_TIN.out
      }
      //SkipLowReadBam(bams_in)
      //disabled for testing
      //LCExtrap( SkipLowReadBam.out.bam )
      LCExtrap( SkipLowReadBam.out.bam.map { sample_id, rg_id, bam, bai -> [sample_id, bam, bai] } )
      
    emit:
      rseqc_tin_logs = rseqc_tin_logs
      rseqc_logs = RSeQC.out
      preseq_lce_logs = LCExtrap.out
      qc_bams = SkipLowReadBam.out.bam
}
