# RNASeq-NF

RNASeq-NF is an NGS analysis pipeline for RNA expression quantification and germline variant calling (GATK4).

The pipeline performs the following tasks.

* Sequence trimming (*TrimGalore*)
* rRNA removal (*SortMeRNA*)
* Mapping and read-group annotation (*STAR*)
* Alignment-QC (*RSeQC, Preseq*)
* PCR duplicate detection (*Sambamba MarkDup*)
* Gene-expression quantification (*HTSeq-count, featureCounts*)
* Transcript quantification (*Salmon*)
* Variant calling (*GATK4*)
* QC report (*MultiQC*)

This implementation is a work in progress and aims to reach feature parity with the [UMCU RNASeq pipeline](https://github.com/UMCUGenetics/RNASeq) while also introducing new features and methods according to developments in the field. Several components have been adapted from the [nf-core rnaseq](https://github.com/nf-core/rnaseq) community pipeline and rewritten in [DSL2](https://www.nextflow.io/docs/edge/dsl2.html) syntax to enable a more modular setup.

## Getting started

### 1. Pipeline setup
# Nextflow
Download the [Nextflow](https://www.nextflow.io/) binary.

# Singularity
Install [Singulariy](https://sylabs.io/guides/3.5/admin-guide/) on the host system. For UMCU users, please follow the instructions on the [HPC wiki](https://wiki.bioinformatics.umcutrecht.nl/bin/view/HPC/SlurmScheduler) on how to use Slurm & Singularity.  

Start an interactive Slurm session on the HPC cluster.
```
srun -n 2 --mem 5G --time 12:00:00 --gres=tmpspace:10G --pty bash
```
The nextflow process needs to run until the analysis (see 4) is finished and all jobs have been scheduled. It is therefore wise to execute the above command within a terminal multiplexer, such as screen or Tmux. Alternatively, the nextflow 

Singulariy environment

Though this should be done by default, ensure that the singularity environment variables point to your $TMPDIR location. In the srun command abouve, we passed 10G to our session.

```
SINGULARITY_LOCALCACHEDIR=${TMPDIR}
SINGULARITY_TMPDIR=${TMPDIR}
```

### 2. Get RNASeq-NF

Clone this repository and submodules. Check out the master branch.

```
git clone --recursive https://github.com/UMCUGenetics/RNASeq-NF.git
```

Ensure that the NextflowModules git branch points to dev-ubec.

### 3. Configuration
3.1 [Resource files](./docs/resources.md) \
3.2 [Settings overview](./docs/settings.md) \
3.3 [Analysis](./docs/config.md) 

## 4. Run RNASeq analysis

Run the pipeline with default setting and genome config 
```
./nextflow run ./RNASeq-NF/main.nf NXF_VER=19.10.0 --fastq_path <fastq_dir> --out_dir <output_dir> --genome_config <path/to/genome.config -profile <slurm,sge>
```

Run the pipeline with custom settings file (`-c`) option.
```
./nextflow run ./RNASeq-NF/main.nf NXF_VER=19.10.0 --fastq_path <fastq_dir> --out_dir <output_dir> --genome_config <path/to/genome.config -c </path/to/myrun.config> -profile <slurm,sge>
```

For local execution, simply omit the -profile parameter.










