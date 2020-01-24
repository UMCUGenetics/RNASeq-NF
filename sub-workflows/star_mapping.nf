include AlignReads from '../../NextflowModules/STAR/2.4.2a/AlignReads.nf' params(params)
include Count from '../../NextflowModules/HTSeq/0.6.0/Count.nf' params(params)
include Index from '../../NextflowModules/Sambamba/0.6.8/Index.nf' params(params)

workflow star_mapping {
    get:
      fastqs
      genome_index
    main:
      /* Run mapping on a per sample per lane basis */
      AlignReads(fastqs,genome_index)
      Index(AlignReads.out.map { sample_id, bams, unmapped, log1, log2, tab -> [sample_id, bams] })
    emit:
      bams = AlignReads.out
      bais = Index.out
}
