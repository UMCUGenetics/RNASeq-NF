# RNASeq-NF

RNASeq-NF is an NGS analysis pipeline for RNA expression quantification and germline variant calling (GATK4).

The pipeline performs the following tasks.

* Read QC (*FastQC*)
* Sequence trimming (*TrimGalore*)
* rRNA removal (*SortMeRNA*)
* Mapping and read-group annotation (*STAR*)
* Alignment-QC (*RSeQC, Preseq*)
* PCR duplicate detection (*Sambamba MarkDup*)
* Gene-expression/biotype quantification ( *featureCounts*)
* Gene-expression normalization ( *edgeR*, *DESeq2*)
* Transcript quantification (*Salmon*)
* Variant calling (*GATK4*)
* QC report (*MultiQC*, *CustomQC*)

This implementation is a work in progress and aims to reach feature parity with the [UMCU RNASeq pipeline](https://github.com/UMCUGenetics/RNASeq) while also introducing new features and methods according to developments in the field. 

Several components have been adapted from the [nf-core rnaseq](https://github.com/nf-core/rnaseq) community pipeline and rewritten in [DSL2](https://www.nextflow.io/docs/edge/dsl2.html) syntax to enable a more modular setup. 

Please refer to the nf-core rnaseq [manual](https://nf-co.re/rnaseq/docs/output) for a description on how to interpret the different output files.

## Getting started

### 1. Pipeline setup
### Nextflow
Download the [Nextflow](https://www.nextflow.io/) binary.

### Singularity
Install [Singularity](https://sylabs.io/guides/3.5/admin-guide/) on the host system. 

For UMCU users, please follow the instructions on the [HPC wiki](https://wiki.bioinformatics.umcutrecht.nl/bin/view/HPC/SlurmScheduler) on how to use Slurm & Singularity.  

Start an interactive Slurm session on the HPC cluster.
```
srun -n 2 --mem 5G --time 12:00:00 --gres=tmpspace:10G --pty bash
```
The nextflow process needs to remain active until the analysis (see 4) is finished and all jobs have been scheduled. It is therefore wise to execute the above command within a terminal multiplexer, such as **screen** or **Tmux**. Alternatively, the command can be embedded within a **sbatch** script. Though this should be done by default, ensure that the singularity environment variables point to $TMPDIR. In the srun command above, we passed 10G to our session.

```
SINGULARITY_LOCALCACHEDIR=${TMPDIR}
SINGULARITY_TMPDIR=${TMPDIR}
```

### 2. Get RNASeq-NF

Clone this repository and ensure that the master branch is checked-out.

```
git clone --recurse-submodules https://github.com/UMCUGenetics/RNASeq-NF.git
```

### 3. Configuration
3.1 [Resource files](./docs/resources.md) \
3.2 [Settings overview](./docs/settings.md) \
3.3 [Analysis](./docs/config.md) 

## 4. Run RNASeq analysis

Before starting the pipeline, ensure that the input fastq files follow the [Illumina Naming Convention](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.html).

Run the pipeline.
```
./nextflow run ./RNASeq-NF/main.nf -c </path/to/run.config> --fastq_path <fastq_dir> --out_dir <output_dir> -profile slurm -resume 
```
For local execution (without HPC backend), simply omit the -profile parameter.










