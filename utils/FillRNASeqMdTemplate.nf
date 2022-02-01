process FillMdTemplate {
    tag {"Fill Markdown Template"}
    label 'FillMdTemplate_1_0'

    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        val(run_name)
        path(md_file)
        path(parameters_file)        
     
    output:        
        path("${run_name}_RNASeq_report.md", emit: md_file)
   
    script:
        """
        # make copy of values file
        tmp_file=${parameters_file}.tmp
        cp ${parameters_file} \${tmp_file}

        # add other parameters from environment 
        echo "value_pipeline=${workflow.manifest.version}" >> \${tmp_file}
        echo "value_nextflow=${workflow.manifest.nextflowVersion}" >> \${tmp_file}
        echo "value_genome=${params.genome}" >> \${tmp_file}
        echo "value_gtf=${params.genome_gtf.split('/')[-1]}" >> \${tmp_file}
        echo "value_fasta=${params.genome_fasta.split('/')[-1]}" >> \${tmp_file}
        if [[ "${params.single_end}" == "true" ]]; then
            echo "value_mode=Single-end" >> \${tmp_file}
        else
            echo "value_mode=Paired-end" >> \${tmp_file}
        fi

        # echo "value_fasta=${params.genome_fasta.split('/')[-1]}" >> \${tmp_file}

        # construct sed command in order to replace mapped key-value pairs
        sed_cmd="sed "
        while read v; do
            if [[ "\$v" != "" ]]; then  
                key=`echo \$v | cut -f1 -d'='`
                value=`echo \$v | cut -f2 -d'='`
                sed_cmd="\${sed_cmd} -e \"s/\${key}/\${value}/\""
            fi
        done < \${tmp_file}

        eval \${sed_cmd} ${md_file} >  ${run_name}_RNASeq_report.md
        """
}

process CreateTmpValueMap {
    tag {"Create_Tmp_Value_Map"}
    label 'CreateTmpValueMap_1_0'

    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        val(run_name)
     
    output:        
        path("${run_name}_valuemap.txt", emit: value_map_file)
   
    script:
        """
        touch ${run_name}_valuemap.txt
        """    
}