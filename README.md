# RNASeq-NF

RNASeq-NF is an NGS analysis pipeline that performs expression quantification on RNASeq data.

The pipeline performs the following tasks.

* Read quality and adapter trimming (TrimGalore)
* ReadQC (FastQC)
* Mapping and read-group annotation (STAR)
* post-QC(RSeQC, Preseq)
* PCR duplicate detection (MarkDup)
* Expression quantification (HTSeq-count)

This implementation is a work in progress and aims to reach feature parity with the [UMCU RNASeq pipeline](https://github.com/UMCUGenetics/RNASeq) while also introducing new features and methods according to development in the field. Several components have been adapted from the [nf-core rnaseq][https://github.com/nf-core/rnaseq] Nextflow community pipeline and rewritten in [DSL2](https://www.nextflow.io/docs/edge/dsl2.html) syntax to enable a more modular setup.





```
module load Java/1.8.0_60
nextflow run RNAseq-NF/main.nf -c RNAseq-NF/test/test-run.config -profile sge
```

