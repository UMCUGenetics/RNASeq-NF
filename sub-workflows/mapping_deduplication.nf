include Markdup from '../NextflowModules/Sambamba/0.7.0/Markdup.nf' params(mem:params.sambambamarkdup.mem, optional:params.markdup.toolOptions)
include AlignReads from '../NextflowModules/STAR/2.7.3a/AlignReads.nf' params(singleEnd:params.singleEnd, optional:params.star.toolOptions)   
include Index from '../NextflowModules/Sambamba/0.7.0/Index.nf' params(params)

workflow markdup_mapping {
    take:
      fastq_files
      star_index
      genome_gtf
    main:
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
