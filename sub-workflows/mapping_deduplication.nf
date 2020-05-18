include Markdup from '../NextflowModules/Sambamba/0.7.0/Markdup.nf' params(mem:params.sambambamarkdup.mem, optional:params.markdup.toolOptions)
   
  

workflow markdup_mapping {
    take:
      bams_in
    main:
      /* Run mapping on a per sample per lane basis */
      Markdup(bams_in)
    emit:
      bams = Markdup.out
}
