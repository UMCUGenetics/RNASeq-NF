process mergeFeatureCounts {
    tag "mergefeaturecounts ${run_id}" 
    label 'mergefeaturecounts'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
    val run_id
    file(count_tables)

    output:
    file("${run_id}_featureCounts_raw.txt")

    shell:
    """
    //To be Implemented
    """
}
