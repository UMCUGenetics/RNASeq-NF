process getExonLenghts {
    tag "getExonLengths ${genome_gtf.baseName}" 
    label 'getExonLenghts'
    shell = ['/bin/bash', '-euo', 'pipefail']
    container = 'quay.io/biocontainers/bioconductor-genomicfeatures:1.38.0--r36_0'

    input:
    file(genome_gtf)

    output:
    file("${genome_gtf.baseName}_exon_gene_sizes.txt")

    script:
    """
    exon_lengths.R ${genome_gtf} ${genome_gtf.baseName}_exon_gene_sizes.txt
    """

}
