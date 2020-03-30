# RNASeq-NF

RNASeq-NF is an NGS analysis pipeline for RNA expression quantification and germline variant calling (GATK4).

The pipeline performs the following tasks.

* Read quality and adapter trimming (*fastp*)
* Mapping and read-group annotation (*STAR*)
* Alignment-QC (*RSeQC, Preseq*)
* PCR duplicate detection (*Sambamba MarkDup*)
* Expression quantification (*HTSeq-count, featureCounts*)
* Variant calling (*GATK4*)

This implementation is a work in progress and aims to reach feature parity with the [UMCU RNASeq pipeline](https://github.com/UMCUGenetics/RNASeq) while also introducing new features and methods according to developments in the field. Several components have been adapted from the [nf-core rnaseq](https://github.com/nf-core/rnaseq) community pipeline and rewritten in [DSL2](https://www.nextflow.io/docs/edge/dsl2.html) syntax to enable a more modular setup.

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

2.1 Parameters

*** Pipeline settings 
```
  --singleEnd True/False if single-end sequencing
  --unstranded True/False if library is unstranded
  --stranded True/False if library is forward-stranded 
  --revstranded True/False
  --hts_count_type htseq-count annotation feature for expression quantification (e.x gene_id)
  --fc_group_features subread featureCount group (e.x gene_id)
  --fc_count_type subread featureCount type (e.x CDS, five_prime_UTR etc.)
  --skipMergeLanes True/False skip merging for fastq files from multiple lanes (Salmon)
```
*** Genome settings ***
```
  --genome_fasta  /path/to/reference.fasta
  --genome_dict /path/to/reference.dict
  --genome_index  /path/to/reference.fai
  --genome_gtf /path/to/annotation.gtf
  --genome_bed12 /path/to/annotation.bed12 (required for RSeQC)
  --genome_known_sites /path/to/snp_sites.vcf (list)
  --scatter_interval_list /path/to/scatter.intervals (GATK4)
```
*** Transcriptome settings **
```
  --star_index /path/to/star index
  --salmon_index /path/to/salmon index
  --gene_len /path/to/feature lengths (required for HTSeq RPKM normalization)
```





2.2 Run Analysis.

```
./nextflow run ./RNASeq-NF/main.nf -c ./RNASeq-NF/conf/your_config.conf --fastq_path <fastq_dir>  --out_dir<output_dir> -profile slurm
```
For local execution, simply omit the -profile parameter. You can provide all necessary arguments either directly to nextflow or setup a config file for convenience (See RNASeq-NF/conf/test-run.config) for an example. 







Replace test-run-config with your own configuration.
```
module load Java/1.8.0_60
nextflow run RNAseq-NF/main.nf -c RNAseq-NF/test/test-run.config -profile {sge,slurm}
```
For local execution, simply remove the -profile argument.


