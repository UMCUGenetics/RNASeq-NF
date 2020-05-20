#!/usr/bin/env nextflow
nextflow.preview.dsl=2
include extractAllFastqFromDir from './NextflowModules/Utils/fastq.nf' params(params)
//Workflows
include pre_processing from './sub-workflows/pre_processing.nf' params(params)
include post_mapping_QC from './sub-workflows/post_mapping_QC.nf' params(params)
include markdup_mapping from './sub-workflows/mapping_deduplication.nf' params(params)
include alignment_free_quant from './sub-workflows/alignment_free_quant.nf' params(params)
include alignment_based_quant from './sub-workflows/alignment_based_quant.nf' params(params)
include multiqc_report from './sub-workflows/multiqc_report.nf' params(params)
include gatk_germline_snp_indel from './sub-workflows/gatk_germline_snp_indel.nf' params(params)
//End workflows                                                                  

//Check minimal resource parameters
if (!params.out_dir) {
   exit 1, "Output directory not found. Please provide the correct path!"
}
if (!params.fastq_path) {
  exit 1, "fastq files not found. Please provide the correct path!"
}
if (!params.genome_fasta) {
  exit 1, "Genome fasta not found. Please provide the correct path!"
} 
if (!params.genome_gtf) {
  exit 1, "Genome GTF not found. Please provide the correct path!"
} 
if (!params.transcripts_fasta && params.runSalmon) {
  exit 1, "Transcript fasta not found. Please provide the correct path!"
} 
if (!params.genome_dict && params.runGATK4_HC) {
    exit 1, "Genome dictionary not found. Please provide the correct path!"
} 
if (!params.genome_index && params.runGATK4_HC) {
    exit 1, "Genome index not found. Please provide the correct path!"
} 
//Start workflow
workflow {
  main :
    //Set run and retrieve input fastqs  
    run_name = params.fastq_path.split('/')[-1]
    fastq_files = extractAllFastqFromDir(params.fastq_path).map { [it[0],it[1],it[4]]}
   //Pipeline log info
    params.version = "Beta"
    log.info """=======================================================
    RNASeq-NF ${params.version}"
    ======================================================="""
    def summary = [:]
    summary['Pipeline Name']  = 'RNASeq-NF'
    summary['Pipeline Version'] = params.version
    summary['Run Name']     = run_name
    summary['Fastq dir']   = params.fastq_path
    summary['Genome config']   = params.genome_config
    summary['Mode']   = params.singleEnd ? 'Single-end' : 'Paired-end'
    summary['Output dir']   = params.out_dir
    summary['Working dir']  = workflow.workDir
    summary['Container Engine'] = workflow.containerEngine
    summary['Current home']   = "$HOME"
    summary['Current user']   = "$USER"
    summary['Current path']   = "$PWD"
    summary['Script dir']     = workflow.projectDir
    summary['Config Profile'] = workflow.profile
    log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
    log.info "========================================="
    // # 1) Pre-processing / QC
    pre_processing ( fastq_files )
    trim_logs = pre_processing.out.trim_logs
    fastqc_logs = pre_processing.out.fastqc_logs
    sortmerna_logs = Channel.empty()
    // Determine final fastqs files
    final_fastqs = pre_processing.out.processed_fastqs 
    //Transform output channels
    if (params.singleEnd) {
        fastqs_transformed = final_fastqs
             .groupTuple(by:0)
             .map { sample_id, rg_ids, reads -> [sample_id, rg_ids[0], reads.flatten().toSorted(), []] }
    } else {
         fastqs_transformed = final_fastqs
              .map{ sample_id, rg_ids, reads -> [sample_id, rg_ids, reads[0], reads[1]] }
              .groupTuple(by:0)
              .map{ sample_id, rg_ids, r1, r2 -> [sample_id, rg_ids[0], r1.toSorted(), r2.toSorted()] }
    }
    // # 2) STAR alignment | Sambamba markdup
    if ( params.runMapping ) {
        mapped = markdup_mapping( fastqs_transformed )
        star_logs = mapped.logs
    } else {
        star_logs = Channel.empty()
    }
    // # 3) Post-mapping QC
    if ( params.runPostQC) {
      if (params.runMapping) {
          post_mapping_QC( mapped.bam_sorted.map { sample_id, rg_id, bam, bai -> [sample_id, bam, bai] } )
          post_qc_logs = post_mapping_QC.out[1].map { it[1] }.mix(post_mapping_QC.out[0].map { it[1] })
      } else {
          exit 1, "PostQC requires alignment step. Please enable runMapping!"
      } 
    } else {
        post_qc_logs = Channel.empty()
    } 
    // # 4) featureCounts 
    if ( params.runFeatureCounts) {
      if (params.runMapping) {
         alignment_based_quant ( run_name, mapped.bam_sorted.map { it[2] } )
         fc_logs = alignment_based_quant.out.fc_summary
        } else {
          exit 1, "featureCounts requires alignment step. Please enable runMapping!"
      } 
    } else {
          fc_logs = Channel.empty()
    }
    // # 5) Salmon    
    if ( params.runSalmon ) {
      alignment_free_quant( fastqs_transformed, run_name )
      salmon_logs = alignment_free_quant.out.logs
    } else {
        salmon_logs = Channel.empty()
    }
    // # 6) GATK4 germline variant calling with optional BQSR
    if ( params.runGATK4_HC ) {
      if ( params.runMapping ) {
          gatk_germline_snp_indel( run_name, mapped.bam_dedup.map {sample_id, rg_id, bam, bai -> 
                                                                 [sample_id, bam, bai] }) 
      } else {
            exit 1, "GATK4 requires alignment, markdup step. Please enable runMapping and runMarkDup!"
      }        
    }
    // # 7) MultiQC report  
    if ( params.runMultiQC ) {
        multiqc_report( run_name,
                        fastqc_logs,
                        trim_logs,
                        sortmerna_logs,
                        star_logs,
                        post_qc_logs,
                        fc_logs,
                        salmon_logs ) 		      
      }

}
