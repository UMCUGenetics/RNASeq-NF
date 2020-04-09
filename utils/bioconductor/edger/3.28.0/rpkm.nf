process rpkm {
    tag "rpkm ${run_id}"
    label 'biconductor_3_20_7'
    label 'biconductor_3_20_7_edger_rpkm'
    
    container = 'quay.io/biocontainers/bioconductor-edger:3.20.7--r3.4.1_0'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
    val run_id
    file(counts)
    val feature_lengths

    output:
    file("${run_id}_${params.tool}_readCounts_RPKM.txt")

    script:
    """
    edgerRpkm.R ${run_id}_${params.tool} ${counts} ${feature_lengths}  
    """

}
