include TrimGalore from '../NextflowModules/TrimGalore/0.6.1/TrimGalore.nf' params(params)

workflow pre_mapping_QC {
    get:
      fastqs
    main:
      /* Run mapping on a per sample per lane basis */
      TrimGalore(fastqs)
    emit:
      TrimGalore.out
      
}
