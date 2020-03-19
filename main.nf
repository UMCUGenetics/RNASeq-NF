#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include './NextflowModules/Utils/fastq.nf' params(params)
include post_mapping_QC from './sub-workflows/post_mapping_QC.nf' params(params)
include markdup_mapping from './sub-workflows/mapping_deduplication.nf' params(params)
include multiqc_report from './sub-workflows/multiqc_report.nf' params(params)
include gatk4_bqsr from './sub-workflows/gatk4_bqsr.nf' params(params)
include gatk4_hc from './sub-workflows/gatk4_hc.nf' params(params)
include SplitIntervals from './NextflowModules/GATK/4.1.3.0/SplitIntervals.nf' params(params)
include SplitNCigarReads from './NextflowModules/GATK/4.1.3.0/SplitNCigarReads.nf' params(params)
include Count from './NextflowModules/HTSeq/0.6.0/Count.nf' params(params)
include AlignReads from './NextflowModules/STAR/2.6.0c/AlignReads.nf' params(params)
include Index from './NextflowModules/Sambamba/0.6.8/Index.nf' params(params)
include gatk4_rnaseq from './sub-workflows/gatk4_rnaseq.nf' params(params)
include Quant from './NextflowModules/Salmon/0.13.1/Quant.nf' params(params)
include Fastp from './NextflowModules/fastp/0.14.1/Fastp.nf' params(params)
include mergeFastqLanes from './NextflowModules/Utils/mergeFastqLanes.nf' params(params)
include mergeHtseqCounts from './NextflowModules/Utils/mergeHtseqCounts.nf' params(params)
include edgerRpkm from './NextflowModules/Utils/edgerRpkm.nf' params(params)

if (!params.fastq_path) {
   exit 1, "fastq directory does not exist. Please provide correct path!"
}
if (!params.out_dir) {
   exit 1, "Output directory not found. Please provide the correct path!"
}

workflow {
  main :
    run_id = "Test_run"
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
//    if (!params.skipMapping && !params.skipMarkDup && !params.skipGATK4) {
//      genome_fasta = Channel
//            .fromPath(params.genome_fasta, checkIfExists: true)
//            .ifEmpty { exit 1, "Genome fasta file not found: ${params.genome_fasta}"}
//      genome_idx = Channel
//            .fromPath(params.genome_index, checkIfExists: true)
//            .ifEmpty { exit 1, "Genome index not found: ${params.genome_index}"}
//      genome_dict = Channel
//            .fromPath(params.genome_dict, checkIfExists: true)
//            .ifEmpty { exit 1, "Genome dictionary not found: ${params.genome_dict}"}
//    }
    if (!params.skipSalmon) {
      salmon_index = Channel
            .fromPath(params.salmon_index, checkIfExists: true)
            .ifEmpty { exit 1, "Transcripts fasta not found: ${params.salmon_index}"}
    }
    if (params.singleEnd) {
      if (!params.skipFastp) {
	    final_fastqs = Fastp(fastq_files)
            .groupTuple(by:0)
            .map { sample_id, rg_ids, json, reads -> [sample_id, rg_ids[0], reads.toSorted(), [], json] }
            
 	  } else {
        final_fastqs = fastq_files
            .groupTuple(by:0)
            .map { sample_id, rg_ids, reads -> [sample_id, rg_ids[0], reads.flatten().toSorted(), []] }
      }
    //Paired-end mode         
    } else {
        if (!params.skipFastp) {
          final_fastqs =  Fastp(fastq_files)
            .map{ sample_id, rg_ids, json, reads -> [sample_id, rg_ids, reads[0], reads[1], json] }
            .groupTuple(by:0)
            .map{ sample_id, rg_ids, r1, r2, json -> [sample_id, rg_ids[0], r1.toSorted(), r2.toSorted(), json] }
        } else {
          final_fastqs = fastq_files
            .map{ sample_id, rg_ids, reads -> [sample_id, rg_ids, reads[0], reads[1]] }
            .groupTuple(by:0)
            .map{ sample_id, rg_ids, r1, r2 -> [sample_id, rg_ids[0], r1.toSorted(), r2.toSorted()] }
        }
    } 
    if (!params.skipMapping) {
      AlignReads(final_fastqs.map { [it[0], it[1], it[2], it[3]] }, genome_index.collect())
      Index(AlignReads.out.map { sample_id, bams, unmapped, log1, log2, tab -> [sample_id, bams] })
      mapped = AlignReads.out.join(Index.out)
    }
    if (!params.skipPostQC && !params.skipMapping) {
      post_mapping_QC(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, bams, bai] }, genome_bed.collect())
    }
    if (!params.skipCount && !params.skipMapping) {
      Count(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, bams, bai] }, genome_gtf.collect())
      //Merge HTSeq counts
      mergeHtseqCounts( "Test_run", Count.out.map { it[1] }.collect())
      //RPKM counts
      edgerRpkm("Test_run", mergeHtseqCounts.out, params.gene_len)

    }
    if (!params.skipMarkDup && !params.skipMapping) {
      markdup_mapping(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, sample_id, bams, bai] })
    }
    if (!params.skipSalmon ) {
      if (!params.skipMergeLanes) {
        Quant ( mergeFastqLanes (final_fastqs), salmon_index.collect() )
      } else if (!params.singleEnd && params.skipMergeLanes) {
	  final_fastqs.view()
          Quant ( final_fastqs.map {sample_id, rg_id, r1, r2, json -> [sample_id, [r1,r2].flatten()] }, salmon_index.collect() )
      } else {
          Quant ( final_fastqs.map {sample_id, rg_id, reads, json -> [sample_id, reads] }, salmon_index.collect() ) 
      }
    } 
    if (!params.skipMapping && !params.skipMarkDup && !params.skipGATK4_HC) {
          //Split ncigar's
          SplitNCigarReads(markdup_mapping.out)
          SplitIntervals( 'no-break', Channel.fromPath( params.scatter_interval_list))
          if (!params.skipGATK4_BQSR) {
              //Perform BSQR
             gatk4_bqsr(SplitNCigarReads.out, SplitIntervals.out.flatten())
             gatk4_hc(gatk4_bqsr.out[0], SplitIntervals.out.flatten())
          } else {
             gatk4_hc(SplitNCigarReads.out, SplitIntervals.out.flatten())
          }
            
    }
    if (!params.skipMultiQC) {
      multiqc_report( final_fastqs.map { it[-1] }, 
		      AlignReads.out.map{ [it[3], it[4]] }, 
                      post_mapping_QC.out[1].map { it[1] }.mix(post_mapping_QC.out[0].map { it[1] }),  
                      Count.out.map { it[1] },
		      gatk4_bqsr.out[1].map {it[1]})
   }

}
