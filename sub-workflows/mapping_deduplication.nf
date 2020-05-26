include Markdup from '../NextflowModules/Sambamba/0.7.0/Markdup.nf' params( mem:params.sambambamarkdup.mem )
include AlignReads from '../NextflowModules/STAR/2.7.3a/AlignReads.nf' params( singleEnd:params.singleEnd, optional:params.options.STAR )   
include Index from '../NextflowModules/Sambamba/0.7.0/Index.nf' params( params )
include GenomeGenerate from '../NextflowModules/STAR/2.7.3a/GenomeGenerate.nf' params( params )

workflow markdup_mapping {
    take:
      fastq_files
      genome_fasta
      genome_gtf

    main:
     //Create STAR index if not present
      if (params.star_index) {
          star_index = Channel
              .fromPath(params.star_index, checkIfExists: true)
              .ifEmpty { exit 1, "STAR index not found: ${params.star_index}"}
      } else if (!params.star_index && params.runMapping) {
          //Generate Genome
          GenomeGenerate ( genome_fasta, genome_gtf )
          star_index = GenomeGenerate.out.star_index
      } 
      /* Run mapping on a per sample per lane basis */
      AlignReads( fastq_files, star_index.collect(), genome_gtf.collect() )
      Index(AlignReads.out.bam_file.map {sample_id, rg_id, bam ->
                                         [sample_id, bam] })
      Markdup(AlignReads.out.bam_file.join(Index.out))
      
    emit:
      bam_sorted = AlignReads.out.bam_file.join(Index.out)
      logs = AlignReads.out.log.mix(AlignReads.out.final_log)
      bam_dedup = Markdup.out
}
