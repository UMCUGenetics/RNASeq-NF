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

git clone --recursive https://github.com/UMCUGenetics/RNASeq-NF.git
```


Make sure you have access to the Singularity containers on the HPC. 
```
/hpc/local/CentOS7/cog_bioinf/nextflow_containers
```
Technically, this should also work with Biocontainers.

***Usage instructions***

Replace test-run-config with your own configuration.
```
module load Java/1.8.0_60
nextflow run RNAseq-NF/main.nf -c RNAseq-NF/test/test-run.config -profile {sge,slurm}
```
For local execution, simply remove the -profile argument.


