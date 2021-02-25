include RSeQC from '../NextflowModules/RSeQC/3.0.1/RSeQC.nf' params( singleEnd:params.singleEnd)
include RSeQC_TIN from '../NextflowModules/RSeQC/3.0.1/RSeQC.nf' params( params )
include LCExtrap from '../NextflowModules/Preseq/2.0.3/LCExtrap.nf' params( optional:params.options.Preseq )
include GtfToGenePred from '../NextflowModules/UCSC/377/GtfToGenePred.nf' params( params )
include GenePredToBed from '../NextflowModules/UCSC/377/GenePredToBed.nf' params( params )

workflow post_mapping_QC {
    take:
      bams_in
      genome_gtf

    main:
        rseqc_logs = Channel.empty()
        rseqc_tin_logs = Channel.empty()
        preseq_lce_logs = Channel.empty()
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
        RSeQC(bams_in, genome_bed.collect())
      
      if (params.runRSeQC_TIN) {
          RSeQC_TIN(bams_in, genome_bed.collect()) 
            rseqc_tin_logs = RSeQC_TIN.out
      }
      LCExtrap(bams_in)
      
    emit:
      rseqc_tin_logs = rseqc_tin_logs
      rseqc_logs = RSeQC.out
      preseq_lce_logs = LCExtrap.out
}
