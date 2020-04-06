### 2.2 Configuration files

**Analysis configuration**

Default settings are stored in `./conf/base.config` and will be automatically set upon execution the pipeline. These defaults can be overwriten by either appending the setting of interest directly to the command-line or storing them in a seperate configuration file (see `./conf/tet_run.config`) and appending them with `-c my_run.config`. 

For example:

*  `./nextflow run ./RNASeq-NF/main.nf -c ./RNASeq-NF/conf/my-run.config ....` (configuration file)
*  `./nextflow run ./RNASeq-NF/main.nf --skipCount=true.....`(without configuration file)

**Genome configuration**

Custom genome settings should be included as a separate configuration file via the the `--genome_config` parameter, either directly or within a configuration file as described above.  For example;

```
includeConfig '../process.config

params {
  star_index = '/hpc/cog_bioinf/GENOMES/STAR/Homo_sapiens.GRCh37'
  transcripts_fasta = '/hpc/cog_bioinf/GENOMES/RSEM/GRCh37/GRCh37.transcripts.fa'
  genome_gtf = '/hpc/cog_bioinf/GENOMES/STAR/Homo_sapiens.GRCh37/Homo_sapiens.GRCh37.74.gtf'
  genome_bed = '/hpc/cog_bioinf/ubec/tools/RSeQC/Homo_sapiens.GRCh37.74.bed12'
  genome_fasta = '/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fa'
  salmon_index = '/hpc/cog_bioinf/GENOMES/RSEM/GRCh37/GRCh37_transcripts_salmon'
  genome_known_sites = ['/hpc/cog_bioinf/common_dbs/GATK_bundle/1000G_phase1.indels.b37.vcf',
  '/hpc/cog_bioinf/common_dbs/GATK_bundle/dbsnp_137.b37.vcf',
  '/hpc/cog_bioinf/common_dbs/GATK_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf']
  scatter_interval_list = '/hpc/cog_bioinf/ubec/tools/Homo_sapiens.GRCh37.GATK.illumina.chromosomes.interval_list'
  gene_len = '
```

In it's most basic form, the pipeline requires the following resource parameters.

* `--genome_fasta` - path to genome sequence (.fasta). The index (.fai) should be stored in the same directory.
* `--genome_gtf` - path to genome annotation (.gtf)
* `--transcripts_fasta` - path to transcript sequences (.fasta) when skipSalmon is set to false.

All other files will be created at run-time and stored in the output folder (reference_genome). Since buidling the STAR index requires a significant amount of memory, ensure that sufficient memory is available on the host system and allocated in process.config. 

**Process configuration**

Runtime specific resources (memory, cpu's) can be configured in the process.config file. These settings should be sufficent for most jobs but can be always be altered if equired.

```
 withLabel : HTSeq_0_11_3_Count {
      time = '24h'
      penv = 'threaded'
      cpus = 2
```
**Nextflow configuration**

Nextflow's base configuration (executors, containerization etc.) are stored in nextflow.config. Ensure that nextflow.config and process.config are always findable by the nextflow runner.

