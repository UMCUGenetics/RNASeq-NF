include SplitNCigarReads from '../NextflowModules/GATK/4.1.3.0/SplitNCigarReads.nf' params(params)

workflow gatk4_rnaseq {
    get:
      bam_dedup
      genome_fasta
      genome_index
      genome_dict
    main:
      SplitNCigarReads(bam_dedup, genome_fasta.collect(), genome_index.collect(), genome_dict.collect())
    emit:
      split_bam = SplitNCigarReads.out
}
