**General pipeline settings**
* `--singleEnd` true/false for single-end sequencing (Default: false)
* `--unstranded` true/false for unstranded library prep (Default: false)
* `--stranded` true/false for forward-stranded library prep (Default: false)
* `--revstranded` true/false for reverse-stranded library prep (Default: true)
* `--hts_group_features` htseq-count annotation feature for expression quantification (Default: gene_id)
* `--hts_count_type` htseq-count annotation feature for expression quantification (Default: exon)
* `--fc_group_features` subread featureCounts group (Default: gene_id)
* `--fc_count_type subread` subread featureCounts type (Default: exon)
* `--normalize_counts` enable edgeR RPKM/CPM normalization for featureCounts (Default: true)
* `--gencode` gencode reference (Default: false).

**Pipeline steps**

* `--skipTrimGalore` Skip trimming with TrimGalore (Default: false)
* `--skipSortMeRna` Skip rRNA filtering with SortMeRNA (Default: false)
* `--skipPostQC` Skip post alignment QC (Default: false)
* `--skipMarkDup` Skip Sambamba markdup (Default: false)
* `--skipHTSeqCount` Skip htseq-count expression quantification (Default: false)
* `--skipFeatureCounts` Skip featureCounts expression quantification (Default: false)
* `--skipMapping` Skip STAR alignment (Default: false)
* `--skipSalmon`Skip Salmon transcript quantification (Default: false)
* `--skipMultiQC` Skip MultiQC quality report (Default: false)
* `--skipGATK4_HC` Skip GATK4 variant calling (Default: true)
* `--skipGATK4_BQSR` Skip GATK4 base quality recalibration (Default: true)

**Reference resources**

* `--genome_fasta` /path/to/reference/genome.fasta
* `--genome_index` /path/to/reference/genome.fasta.fai
* `--genome_dict` /path/to/reference/genome.dict
* `--genome_gtf` /path/to/reference/annotation.gtf
* `--transcripts_fasta` /path/to/reference/transcripts.fasta (required for Salmon)
* `--genome_bed` path/to/reference/annotation.bed12 
* `--star_index` path/to/star_index 
* `--salmon_index` path/to/salmon_index  
* `--genome_known_sites` path/to/snp_sites.vcf (optional, GATK4 BQSR) 
* `--scatter_interval_list` path/to/scatter.interval_list (required for GATK4) 
* `--rRNA_database_manifest` path/to/rRNA database file (required for SortMeRNA) 

