# RNASeq-NF

RNASeq-NF is an NGS analysis pipeline for RNA expression quantification and germline variant calling (GATK4).

The pipeline performs the following tasks.

* Read quality and adapter trimming (*fastp*)
* Mapping and read-group annotation (*STAR*)
* Alignment-QC (*RSeQC, Preseq*)
* PCR duplicate detection (*Sambamba MarkDup*)
* Gene-expression quantification (*HTSeq-count, featureCounts*)
* Transcript expression-quantification (*Salmon)
* Variant calling (*GATK4*)
* QC report (*MultiQC*)

This implementation is a work in progress and aims to reach feature parity with the [UMCU RNASeq pipeline](https://github.com/UMCUGenetics/RNASeq) while also introducing new features and methods according to developments in the field. Several components have been adapted from the [nf-core rnaseq](https://github.com/nf-core/rnaseq) community pipeline and rewritten in [DSL2](https://www.nextflow.io/docs/edge/dsl2.html) syntax to enable a more modular setup.

## Documentation

1. [Installation](./docs/installation.md) 
2. [Configuration](./docs/config.md) \
2.1 [Settings](./docs/settings.md) \
2.2 [Resources](./docs/reference.md) 
3. [Run RNASeq-NF](./docs/running.md) 









