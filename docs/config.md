### 2.2 Configuration files

**1. Analysis configuration**

Default [settings](./settings.md) are stored in `./conf/base.config` and will be automatically set upon execution of the pipeline. These defaults can be overwriten by either appending the setting of interest directly as command-line arguments or storing them in a seperate configuration file. An example of such a configuration file can be found in  `./conf/tet_run.config`. 

Runtime specific resources (memory, cpu's) should be sufficent for most jobs but can be always be altered if required. 

```
 withLabel : HTSeq_0_11_3_Count {
      time = '24h'
      penv = 'threaded'
      cpus = 2
```

The new configuration can be appended with the  `-c` option. For example;

*  `./nextflow run ./RNASeq-NF/main.nf -c ./RNASeq-NF/conf/my-run.config ....` (configuration file)
*  `./nextflow run ./RNASeq-NF/main.nf --skipCount=true.....`(without configuration file)

**2. Resource configuration**

Genomic resources files should always be included as a separate configuration file via the the `--genome_config` parameter, either directly or within a configuration file as described above. 

In it's most basic form, the pipeline requires the following resource parameters.

* `--genome_fasta` - path to genome sequence (.fasta). The index (.fai) should be stored in the same directory.
* `--genome_gtf` - path to genome annotation (.gtf)

Transcript quantification with Salmon requires transcript sequences in .fasta format (In case no Salmon index is provided)
* `--transcripts_fasta` - path to transcript sequences (.fasta) when skipSalmon is set to false.


I
GATK requires a sequence dictionary (.dict) and fasta index (.fai) to be present in the same directionary as the reference sequence. Please prepare them before running the pipeline if applicable as described in step 2 and 3 of the [GATK reference preparation guide](https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk)

Known variants in vcf format (optional) should be included in the genome config as shown below.

```
 genome_known_sites = ['/hpc/cog_bioinf/common_dbs/GATK_bundle/1000G_phase1.indels.b37.vcf',
  '/hpc/cog_bioinf/common_dbs/GATK_bundle/dbsnp_137.b37.vcf',
  '/hpc/cog_bioinf/common_dbs/GATK_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf']
```

All remaining files will be created at run-time and stored in the output folder (reference_resources). Since buidling the STAR index requires a significant amount of memory, ensure that sufficient memory is available on the host system and allocated in process.config. 

Example genome configuration files can be found `./resources/test_run`.

**Nextflow configuration**

Nextflow's base configuration (executors, containerization etc.) are stored in the nextflow.config file.

Set the cacheDir to a location with enough space to store the singularity images. 

`cacheDir = '/path/to/my/cache/dir`

Ensure that nextflow.config and process.config are always findable by the nextflow runner.

