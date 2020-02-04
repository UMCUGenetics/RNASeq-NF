include Index from '../NextflowModules/Sambamba/0.6.8/Index.nf' params(params)
include MarkDup from '../NextflowModules/Sambamba/0.6.8/MarkDup.nf' params(params)

workflow markdup_mapping {
    get:
      bams_in
    main:
      /* Run mapping on a per sample per lane basis */
      MarkDup(bams_in)
    emit:
      bams = MarkDup.out
}
