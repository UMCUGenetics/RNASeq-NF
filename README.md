# RNASeq-NF

RNASeq-NF is an NGS analysis pipeline for RNA expression quantification and germline variant calling (GATK4).

The pipeline performs the following tasks.

* Read quality and adapter trimming (*fastp*)
* Mapping and read-group annotation (*STAR*)
* Alignment-QC (*RSeQC, Preseq*)
* PCR duplicate detection (*Sambamba MarkDup*)
* Expression quantification (*HTSeq-count, featureCounts*)
* Variant calling (*GATK4*)
* QC report (*MultiQC*)

This implementation is a work in progress and aims to reach feature parity with the [UMCU RNASeq pipeline](https://github.com/UMCUGenetics/RNASeq) while also introducing new features and methods according to developments in the field. Several components have been adapted from the [nf-core rnaseq](https://github.com/nf-core/rnaseq) community pipeline and rewritten in [DSL2](https://www.nextflow.io/docs/edge/dsl2.html) syntax to enable a more modular setup.

## Remarks ###
* Set skipMultiQC=true. There are currently some issues with MultiQC and optional steps (work in progress) 
* Reference file generation is currently not automated. Look in /resources/prep_gene for an example on how to prepare a subset of the reference files.

## 1. Installation

### 1.1 Nextflow 
Install [Nextflow](https://www.nextflow.io/) 

### 1.2 Singularity
Install [Singulariy](https://sylabs.io/guides/3.5/admin-guide/) on the host system. Required for biocontainers.

### 1.3 Clone this repository and submodules

```
git clone --recursive https://github.com/UMCUGenetics/RNASeq-NF.git
```

## 2.Usage

### 2.1 Parameters

**Pipeline settings** 
* `--singleEnd` True/False for single-end sequencing
* `--unstranded` True/False for unstranded library prep
* `--stranded` True/False for forward-stranded library prep
* `--revstranded` True/False for reverse-stranded library prep
* `--hts_count_type` htseq-count annotation feature for expression quantification (<em> gene_id, transcript_id etc.</em>)
* `--fc_group_features` subread featureCount group (e.x gene_id)
* `--fc_count_type subread` featureCount type (<em> CDS, five_prime_UTR etc.</em>)
* `--skipMergeLanes` Skip merging Fastq files from multiple lanes (Salmon)
* `--skipPostQC` Skip post alignment QC (RSeQC, Preseq)
* `--skipMarkDup` Skip Sambamba markdup
* `--skipCount` Skip HTSeq/featureCounts read quantification.
* `--skipMapping` Skip STAR alignment
* `--skipSalmon`Skip Salmon transcript quantification
* `--skipFastp` Skip trimming with fastp
* `--skipMultiQC` Skip MultiQC quality report
* `--skipGATK4_HC` Skip GATK4 variant calling
* `--skipGATK4_BQSR` Skip GATK4 base quality recalibration
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


### 2.2 Configuration files

**Resource configuration**

You can provide all necessary parameters either directly to nextflow on the command-line or setup config files for your convenience. For example, store all genomic resource settings in a single config file and include it in your run specific configuration along with process.config (see 2.3). 
```
includeConfig '../process.config'
includeConfig '../resources/UMCU_GRCh37.config'
```

You can find example of genomic resource and run specific configuration files in the resources (<em>/resources</em>) and conf (<em>/conf/test-run.config</em>) folders.



**Process configuration**

Runtime specific resources (memory, cpu's) can be configured in process.config. Furthermore, advanced parameters can be set via params.<tool>.toolOptions for individual components.  
```
 withLabel : HTSeq_0_11_3_Count {
      params.count.mem = '25G'
      params.count.toolOptions = '-m union -r pos'
      time = '24h'
      penv = 'threaded'
      cpus = 2
```
  
**Nextflow configuration**

Nextflow's base configuration settings (executors, containerization etc.) are stored in nextflow.config.

## 2.2 Run Analysis.

```
module load Java/1.8.0_60
./nextflow run ./RNASeq-NF/main.nf -c <your_run.config>? --fastq_path <fastq_dir>  --out_dir <output_dir> -profile slurm
```
For local execution, simply omit the -profile parameter. 







