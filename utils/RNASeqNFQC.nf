process RNASeqNFQC {
    tag "RNASeqNFQC ${run_id}"
    label 'rnaseqnfqc_1_0_0'
    label 'rnaseqnfqc_1_0_0'
    
    container = 'file:///hpc/compgen/users/tschafers/container/rb-5f115574487994a95c1139a3_latest.sif'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       val(run_id)
       file(qc_files)

    output:
       file("*.html")
   
    script:
       """
       Rscript -e "rmarkdown::render(\"/bin/RNASeqNF_QC.Rmd\", params = list(input =\"\$PWD\"))"  
       """

}
