process SkipLowReadBam {
    tag {"SkipLowReadBam ${sample_id} "}
    label 'SkipLowReadBam_1_0'
    container = 'quay.io/biocontainers/samtools:1.10--h9402c20_2'
    //disabled because otherwise head for some reason exits with error 141
    //shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(val(sample_id), val(rg_id), path(bam_file), path(bai_file))

    output:
        tuple(val(sample_id), val(rg_id), path(bam_file), path(bai_file), optional:true, emit: bam)

    script:
    """
        #get nr of reads (up to a maximum of params.minReadcount)
        count_value=`samtools view ${bam_file} | head -${params.minReadcount} | wc -l`

        if [[ \$count_value -lt ${params.minReadcount} ]]; then
            mv ${bam_file} ${bam_file}.skip
            mv ${bai_file} ${bai_file}.skip
        fi
    """
}
