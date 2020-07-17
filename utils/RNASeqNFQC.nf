process RNASeqNFQC {
    tag "RNASeqNFQC ${run_id}"
    label 'rnaseqnfqc_1_0_0'
    label 'rnaseqnfqc_1_0_0'
    
    container = 'library://tscha/remote-builds/rb-5f115574487994a95c1139a3'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       val(run_id)
       file(rmd_template)
       file(qc_files)

    output:
       file("*.html")
   
    script:
       """
       Rscript -e "rmarkdown::render(\\"RNASeqNF_QC.Rmd\\", params = list(input =\\"\$PWD\\"))"  
       """

}
