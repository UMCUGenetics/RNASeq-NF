**General pipeline settings**

* `--singleEnd` true/false for single-end sequencing (Default: false)
* `--unstranded` true/false for unstranded library prep (Default: true)
* `--stranded` true/false for forward-stranded library prep (Default: false)
* `--revstranded` true/false for reverse-stranded library prep (Default: false)
* `--hts_group_features` htseq-count annotation feature for expression quantification (Default: gene_id)
* `--hts_count_type` htseq-count annotation feature for expression quantification (Default: exon)
* `--fc_group_features` subread featureCounts group (Default: gene_id)
* `--fc_count_type subread` subread featureCounts type (Default: exon)
* `--fc_group_features_type` GTF biotype field for subread featureCounts (Default: gene_biotype)
* `--normalize_counts` enable edgeR RPKM/CPM normalization for featureCounts (Default: true)
* `--gencode` gencode reference (Default: false).

**Pipeline steps**

* `--runTrimGalore` Read trimming with TrimGalore (Default: true)
* `--runSortMeRna` rRNA filtering with SortMeRNA (Default: true)
* `--runPostQC` Alignment QC with RSeQC,Preseq (Default: true)
* `--runMarkDup` Sambamba markdup (Default: true)
* `--runHTSeqCount` Expression quantification with htseq-count (Default: false)
* `--runFeatureCounts` Expression quantification with featureCounts (Default: true)
* `--runMapping` Read alignment with STAR (Default: true)
* `--runSalmon` Alignment-free transcript quantification with Salmon (Default: true)
* `--runMultiQC` MultiQC report (Default: true)
* `--runGATK4_HC` GATK4 germline variant calling (Default: false)
* `--runGATK4_BQSR` GATK4 base quality score recalibration (Default: false)

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

