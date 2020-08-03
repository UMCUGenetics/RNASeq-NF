process RNASeqNFQC {
    tag "RNASeqNFQC ${run_id}"
    label 'rnaseqnfqc_1_0_0'
    label 'rnaseqnfqc_1_0_0'
    
    container = 'library://tscha/default/rnaseq-nf_qc:v1.0.0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       val(run_id)
       file(qc_files)
       file(rmd_file)
    output:
       file("*.html")
   
    script:
       """
       Rscript -e "rmarkdown::render(\\"RNASeqNF_QC.Rmd\\", params = list(input =\\"\$PWD\\"), output_file =\\"\$PWD/${run_id}_cqc.html\\")"  
       """

}
