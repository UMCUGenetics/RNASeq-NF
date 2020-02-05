# RNASeq-NF

RNASeq-NF is an NGS analysis pipeline that performs expression quantification on RNASeq data.

The pipeline performs the following tasks.

* Read quality and adapter trimming (TrimGalore)
* ReadQC (FastQC)
* Mapping and read-group annotation (STAR)
* post-QC(RSeQC, Preseq)
* PCR duplicate detection (MarkDup)
* Expression quantification (HTSeq-count)

The implementation is work in progress and aims to reach feature parity with the [UMCU RNASeq pipeline](https://github.com/UMCUGenetics/RNASeq).



RNAseq expression quantification pipeline (prototype). Uses DSL2 syntax.


```
module load Java/1.8.0_60
nextflow run RNAseq-NF/main.nf -c RNAseq-NF/test/test-run.config -profile sge
```

