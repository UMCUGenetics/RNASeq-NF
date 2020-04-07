### 2.2 Configuration files

**Analysis configuration**

Default [settings](./settings.md) are stored in `./conf/base.config` and will be automatically set upon execution of the pipeline. These defaults can be overwriten by either appending the setting of interest directly as command-line arguments or storing them in a seperate configuration file. An example of such a configuration file can be found in  `./conf/tet_run.config`. The new configuration can be appended with the  `-c` option. For example;

*  `./nextflow run ./RNASeq-NF/main.nf -c ./RNASeq-NF/conf/my-run.config ....` (configuration file)
*  `./nextflow run ./RNASeq-NF/main.nf --skipCount=true.....`(without configuration file)

**Resource configuration**



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

