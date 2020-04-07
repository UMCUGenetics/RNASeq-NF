### 2.2 Configuration files

**1. Analysis configuration**

Default [settings](./settings.md) are stored in `./conf/base.config` and will be automatically set upon execution of the pipeline. These defaults can be overwriten by either appending the setting of interest directly as command-line arguments or storing them in a seperate configuration file. An example of such a configuration file can be found in  `./conf/tet_run.config`. The new configuration can be appended with the  `-c` option. For example;

*  `./nextflow run ./RNASeq-NF/main.nf -c ./RNASeq-NF/conf/my-run.config ....` (configuration file)
*  `./nextflow run ./RNASeq-NF/main.nf --skipCount=true.....`(without configuration file)

**2. Resource configuration**

Genomic resources files should always be included as a separate configuration file via the the `--genome_config` parameter, either directly or within a configuration file as described above. 

In it's most basic form, the pipeline requires the following resource parameters.

* `--genome_fasta` - path to genome sequence (.fasta). The index (.fai) should be stored in the same directory.
* `--genome_gtf` - path to genome annotation (.gtf)
* `--transcripts_fasta` - path to transcript sequences (.fasta) when skipSalmon is set to false.

All other files will be created at run-time and stored in the output folder (reference_genome). Since buidling the STAR index requires a significant amount of memory, ensure that sufficient memory is available on the host system and allocated in process.config. 

Example genome configuration files can be found `./resources/test_run`.

**Process configuration**

Runtime specific resources (memory, cpu's) can be configured in the process.config file. These settings should be sufficent for most jobs but can be always be altered if equired. The new configuration can be appended with the  `-c` option.

```
 withLabel : HTSeq_0_11_3_Count {
      time = '24h'
      penv = 'threaded'
      cpus = 2
```
**Nextflow configuration**

Nextflow's base configuration (executors, containerization etc.) are stored in nextflow.config. Specify the Singularity chache dir on the host system. All images will be stored at that location.

`cacheDir = '/path/to/my/cache/dir`

Ensure that nextflow.config and process.config are always findable by the nextflow runner.

