include GenomeGenerate from '../NextflowModules/STAR/2.6.0c/GenomeGenerate.nf' params(params)
include Index from '../NextflowModules/Salmon/0.13.1/Index.nf' params( optional:params.salmonindex.toolOptions )
include GtfToGenePred from '../NextflowModules/ucsc/377/gtfToGenePred/GtfToGenePred.nf' params(params)
include GenePredToBed from '../NextflowModules/ucsc/377/genePredToBed/GenePredToBed.nf' params(params)

workflow prep_genome {
    get:
      genome_gtf
      genome_fasta
      transcript_fasta
    main:
      GenomeGenerate( genome_fasta, genome_gtf )
      Index( transcript_fasta )
      GtfToGenePred( genome_gtf )
      GenePredToBed( GtfToGenePred.out )
    emit:
      GenomeGenerate.out
      Index.out
      GenePredToBed.out
}
