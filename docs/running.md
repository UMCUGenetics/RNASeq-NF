
Run the pipeline. For local execution, simply omit the -profile parameter.

```
module load Java/1.8.0_60
export NXF_VER=19.10.0
./nextflow run ./RNASeq-NF/main.nf -c <your_run.config>? --fastq_path <fastq_dir>  --out_dir <output_dir> -profile slurm
```
