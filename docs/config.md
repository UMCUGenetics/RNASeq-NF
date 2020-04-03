### 2.2 Configuration files

**Resource configuration**

All necessary parameters can either be provided directly to nextflow on the command-line or setup config files for your convenience. For example, define all genomic resource settings in a single config f$
```
includeConfig '../process.config'
includeConfig '../resources/UMCU_GRCh37.config'
```

You can find an example of genomic resource and run specific configuration files in the resources (<em>/resources</em>) and conf (<em>/conf/test-run.config</em>) folders.

**Process configuration**

Runtime specific resources (memory, cpu's) can be configured in the process.config file. Furthermore, advanced parameters can be set via params.<tool>.toolOptions for individual processes.
```
 withLabel : HTSeq_0_11_3_Count {
      params.count.mem = '25G'
      params.count.toolOptions = '-m union -r pos'
      time = '24h'
      penv = 'threaded'
      cpus = 2
```

**Nextflow configuration**

Nextflow's base configuration settings (executors, containerization etc.) are stored in nextflow.config. Ensure that nextflow.config and process.config are always findable by the nextflow runner.

