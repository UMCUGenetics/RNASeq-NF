Without configuration file
```
module load Java/1.8.0_60
export NXF_VER=19.10.0
./nextflow run ./RNASeq-NF/main.nf --fastq_path <fastq_dir> --out_dir <output_dir> --genome_config <path/to/genome.config <additional options>
```

With configuration file
```
module load Java/1.8.0_60
export NXF_VER=19.10.0
./nextflow run ./RNASeq-NF/main.nf --fastq_path <fastq_dir> --out_dir <output_dir> --genome_config <path/to/genome.config -c </path/to/myrun.config> <additional options>
```

For local execution, simply omit the -profile parameter.
