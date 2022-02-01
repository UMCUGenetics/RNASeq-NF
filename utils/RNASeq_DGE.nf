process RNASeqNFDGE {
    tag "RNASeqNFDGE ${run_id}"
    label 'rnaseqnfdge_1_0_0'
    label 'rnaseqnfdge_1_0_0'
    
    container = 'library://fmul/default/dge_pdf_nextflow:v1.0.0'

    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
       val(run_id)
       path(count_file)
       path(md_file)
       path(comparisons_file)
       path(dge_template)
    output:
        tuple(path("${run_id}_DGE_Report.pdf"), path("DE"), emit: dge_report)

    script:
        """
        Rscript -e "rmarkdown::render(\\"${dge_template}\\", params = list(comparisons_file =\\"\$PWD/${comparisons_file}\\", md_file =\\"\$PWD/${md_file}\\", count_file =\\"\$PWD/${count_file}\\"), output_file=\\"\$PWD/${run_id}_DGE_Report.pdf\\")"
        """

}
