#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include './NextflowModules/Utils/fastq.nf' params(params)
include star_mapping from './workflows/star_mapping.nf' params(params) 
include pre_mapping_QC from './workflows/pre_mapping_QC.nf' params(params)
include post_mapping_QC from './workflows/post_mapping_QC.nf' params(params)
include markdup_mapping from './workflows/mapping_deduplication.nf' params(params)


workflow {
  main :
    genome_index = Channel.fromPath(params.star_index)
    genome_bed = Channel.fromPath(params.genome_bed)
    fastq_files = extractFastqFromDir(params.fastq_path)
    if (params.singleEnd) {
	trimmed = pre_mapping_QC(fastq_files)
	merged = trimmed.groupTuple(by:0).map { sample_id, rg_ids, reads, logs, fqc -> [sample_id, rg_ids[0], reads, []]}
 	mapped = star_mapping(merged, genome_index.collect())
        mixed = mapped.bams.join(mapped.bais)
        post_mapping_QC(mixed.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, bams, bai] },genome_bed.collect())
        markdup_mapping(mixed.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, sample_id,  bams, bai] })
    }
    

    //trimmed_fq = pre_mapping_QC(fastq_files)
    //             .map{ sample_id, rg_ids, logs, reads -> [sample_id, rg_ids, reads[0], reads[1]] }
    //             .groupTuple(by:0)
    //             .map{ sample_id, rg_ids, r1, r2 -> [sample_id, rg_ids[0], r1, r2] }       
    
    
  

  publish:
    pre_mapping_QC.out to: "${params.out_dir}/PRE-QC/trimmed", mode: 'copy' 
    star_mapping.out to: "${params.out_dir}/mapping/STAR", mode: 'copy'   
    post_mapping_QC.out to: "${params.out_dir}/POST-QC/RSeQC", mode: 'copy'
    markdup_mapping.out to: "${params.out_dir}/mapping/MarkDup", mode: 'copy'
}
