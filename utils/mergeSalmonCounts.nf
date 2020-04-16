process mergeSalmonCounts {
    tag "mergeSalmonCounts ${run_id}" 
    label 'mergeSalmonCounts'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
    val run_id
    file(count_tables)

    output:
    tuple file("${run_id}_Salmon_counts_raw.txt"), file("${run_id}_Salmon_counts_tpm.txt")

    script:
    """
    mergeSalmonCounts.R \$PWD $run_id 
    """

}

