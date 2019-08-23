# RNAseq-NF
RNAseq expression quantification pipeline (prototype)

To run this pipeline, make sure youre using the DSL2 development release (https://github.com/nextflow-io/nextflow/releases/tag/v19.05.0-edge).


```
module load Java/1.8.0_60
sh nextflow run RNAseq-Reloaded/main.nf --outdir /path/to/out/dir --samplesheet /path/to/samplesheet.tsv -profile sge
```

