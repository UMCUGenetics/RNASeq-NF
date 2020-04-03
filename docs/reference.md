**Genome settings**

* `--genome_fasta`  /path/to/reference.fasta
* `--genome_dict` /path/to/reference.dict
* `--genome_index`  /path/to/reference.fai
* `--genome_gtf` /path/to/annotation.gtf
* `--genome_bed12` /path/to/annotation.bed12 (required for RSeQC)
* `--genome_known_sites` /path/to/snp_sites.vcf (list)
* `--scatter_interval_list` /path/to/scatter.intervals (GATK4)

**Transcriptome settings**

* `--star_index` /path/to/star transcriptome index
*  `--salmon_index` /path/to/salmon index
*  `--gene_len` /path/to/gene lengths for RPKM normalization (htseq-count)


