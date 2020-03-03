include SplitIntervals from '../NextflowModules/GATK/4.1.3.0/SplitIntervals.nf' params(params)
include SplitNCigarReads from '../NextflowModules/GATK/4.1.3.0/SplitNCigarReads.nf' params(params)


workflow gatk4_rnaseq {
    get:
      bam_dedup
      genome_fasta
      genome_index
      genome_dict
    main:
      SplitIntervals( 'no-break', Channel.fromPath( params.scatter_interval_list))
      SplitNCigarReads(bam_dedup.combine(SplitIntervals.out.flatten()) , genome_fasta.collect(), genome_index.collect(), genome_dict.collect())
    emit:
      split_bam = SplitNCigarReads.out
}
