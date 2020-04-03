**Pipeline settings**
* `--singleEnd` True/False for single-end sequencing
* `--unstranded` True/False for unstranded library prep
* `--stranded` True/False for forward-stranded library prep
* `--revstranded` True/False for reverse-stranded library prep
* `--hts_count_type` htseq-count annotation feature for expression quantification (<em> gene_id, transcript_id etc.</em>)
* `--fc_group_features` subread featureCount group (e.x gene_id)
* `--fc_count_type subread` featureCount type (<em> CDS, five_prime_UTR etc.</em>)
* `--skipMergeLanes` Skip merging fastq files from multiple lanes (Salmon).
* `--skipPostQC` Skip post alignment QC (RSeQC, Preseq)
* `--skipMarkDup` Skip Sambamba markdup
* `--skipCount` Skip HTSeq/featureCounts expression quantification.
* `--skipMapping` Skip STAR alignment
* `--skipSalmon`Skip Salmon transcript quantification
* `--skipFastp` Skip trimming with fastp
* `--skipMultiQC` Skip MultiQC quality report
* `--skipGATK4_HC` Skip GATK4 variant calling
* `--skipGATK4_BQSR` Skip GATK4 base quality recalibration
