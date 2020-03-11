#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include './NextflowModules/Utils/fastq.nf' params(params)
include post_mapping_QC from './sub-workflows/post_mapping_QC.nf' params(params)
include markdup_mapping from './sub-workflows/mapping_deduplication.nf' params(params)
include Count from './NextflowModules/HTSeq/0.6.0/Count.nf' params(params)
include AlignReads from './NextflowModules/STAR/2.6.0c/AlignReads.nf' params(params)
include Index from './NextflowModules/Sambamba/0.6.8/Index.nf' params(params)
include gatk4_rnaseq from './sub-workflows/gatk4_rnaseq.nf' params(params)
include Quant from './NextflowModules/Salmon/0.13.1/Quant.nf' params(params)
include Fastp from './NextflowModules/fastp/0.14.1/Fastp.nf' params(params)
include mergeFastqLanes from './NextflowModules/Utils/mergeFastqLanes.nf' params(params)

if (!params.fastq_path) {
   exit 1, "fastq directory does not exist. Please provide correct path!"
}
if (!params.out_dir) {
   exit 1, "Output directory not found. Please provide the correct path!"
}

workflow {
  main :
    fastq_files = extractFastqFromDir(params.fastq_path)
    if (!params.skipMapping) {
      genome_index = Channel
            .fromPath(params.star_index, checkIfExists: true)
            .ifEmpty { exit 1, "STAR index not found: ${params.star_index}"}
    }
    if (!params.skipPostQC && !params.skipMapping) {
      genome_bed = Channel
            .fromPath(params.genome_bed, checkIfExists: true)
            .ifEmpty { exit 1, "Bed12 file not found: ${params.genome_bed}"}
    }
    if (!params.skipCount && !params.skipMapping) {
      genome_gtf = Channel
            .fromPath(params.genome_gtf, checkIfExists: true)
            .ifEmpty { exit 1, "GTF file not found: ${params.genome_gtf}"}
    }
    if (!params.skipMapping && !params.skipMarkDup && !params.skipGATK4) {
      genome_fasta = Channel
            .fromPath(params.genome_fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Genome fasta file not found: ${params.genome_fasta}"}
      genome_idx = Channel
            .fromPath(params.genome_index, checkIfExists: true)
            .ifEmpty { exit 1, "Genome index not found: ${params.genome_index}"}
      genome_dict = Channel
            .fromPath(params.genome_dict, checkIfExists: true)
            .ifEmpty { exit 1, "Genome dictionary not found: ${params.genome_dict}"}
    }
    if (!params.skipSalmon) {
      salmon_index = Channel
            .fromPath(params.salmon_index, checkIfExists: true)
            .ifEmpty { exit 1, "Transcripts fasta not found: ${params.salmon_index}"}
    }
    if (params.singleEnd) {
      if (!params.skipFastp) {
	    final_fastqs = Fastp(fastq_files)
            .groupTuple(by:0)
            .map { sample_id, rg_ids, json, reads -> [sample_id, rg_ids[0], reads.toSorted(), []] }
            
 	  } else {
        final_fastqs = fastq_files
            .groupTuple(by:0)
            .map { sample_id, rg_ids, reads -> [sample_id, rg_ids[0], reads.flatten().toSorted(), []] }
      }
    //Paired-end mode         
    } else {
        if (!params.skipFastp) {
          final_fastqs =  Fastp(fastq_files)
            .map{ sample_id, rg_ids, json, reads -> [sample_id, rg_ids, reads[0], reads[1]] }
            .groupTuple(by:0)
            .map{ sample_id, rg_ids, r1, r2 -> [sample_id, rg_ids[0], r1.toSorted(), r2.toSorted()] }
        } else {
          final_fastqs = fastq_files
            .map{ sample_id, rg_ids, reads -> [sample_id, rg_ids, reads[0], reads[1]] }
            .groupTuple(by:0)
            .map{ sample_id, rg_ids, r1, r2 -> [sample_id, rg_ids[0], r1.toSorted(), r2.toSorted()] }
        }
    } 
    if (!params.skipMapping) {
      temp = AlignReads(final_fastqs, genome_index.collect())
      bais = Index(AlignReads.out.map { sample_id, bams, unmapped, log1, log2, tab -> [sample_id, bams] })
      mapped = temp.join(bais)
      //AddOrReplaceReadGroups( AlignReads.out.map { sample_id, rg_id,  bams, unmapped, log1, log2, tab -> [sample_id, rg_id, bams] } )

    }
    if (!params.skipPostQC && !params.skipMapping) {
      post_mapping_QC(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, bams, bai] }, genome_bed.collect())
    }
    if (!params.skipCount && !params.skipMapping) {
      Count(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, bams, bai] }, genome_gtf.collect())
    }
    if (!params.skipMarkDup && !params.skipMapping) {
      markdup_mapping(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, sample_id, bams, bai] })
    }
    if (!params.skipSalmon ) {
      if (!params.skipMergeLanes) {
        Quant ( mergeFastqLanes (final_fastqs), salmon_index.collect() )
      } else if (!params.singleEnd && params.skipMergeLanes) {
          Quant ( final_fastqs.map {sample_id, rg_id, r1, r2 -> [sample_id, [r1,r2].flatten()] }, salmon_index.collect() )
      } else {
          Quant ( final_fastqs.map {sample_id, rg_id, reads -> [sample_id, reads] }, salmon_index.collect() ) 
      }
    } 
    if (!params.skipMapping && !params.skipMarkDup && !params.skipGATK4) {
          gatk4_rnaseq(markdup_mapping.out, genome_fasta, genome_idx, genome_dict)
         
    }
}
