#!/usr/bin/env nextflow

nextflow.preview.dsl=2
include TrimGalore from './NextflowModules/TrimGalore/0.6.5/TrimGalore.nf' params( optional: params.trimgalore.toolOptions, 
                                                                                   singleEnd: params.singleEnd )
include GenomeGenerate from './NextflowModules/STAR/2.7.3a/GenomeGenerate.nf' params(params)
include Index as SalmonIndex from './NextflowModules/Salmon/1.2.1/Index.nf' params( gencode: params.gencode, optional: '' )
include GtfToGenePred from './NextflowModules/UCSC/377/GtfToGenePred/GtfToGenePred.nf' params(params)
include GenePredToBed from './NextflowModules/UCSC/377/GenePredToBed/GenePredToBed.nf' params(params)
include CreateIntervalList from './NextflowModules/Utils/CreateIntervaList.nf' params(params)
include extractAllFastqFromDir from './NextflowModules/Utils/fastq.nf' params(params)
include post_mapping_QC from './sub-workflows/post_mapping_QC.nf' params(params)
include markdup_mapping from './sub-workflows/mapping_deduplication.nf' params(params)
include multiqc_report from './sub-workflows/multiqc_report.nf' params(params)
include SortMeRNA from './NextflowModules/SortMeRNA/4.2.0/SortMeRNA.nf' params(singleEnd:params.singleEnd)
include SplitIntervals from './NextflowModules/GATK/4.1.3.0/SplitIntervals.nf' params(optional: params.splitintervals.toolOptions)
include gatk4_bqsr from './sub-workflows/gatk4_bqsr.nf' params(params)
include gatk4_hc from './sub-workflows/gatk4_hc.nf' params(params)
include SplitNCigarReads from './NextflowModules/GATK/4.1.3.0/SplitNCigarReads.nf' params(genome_fasta:params.genome_fasta)
include Count from './NextflowModules/HTSeq/0.11.3/Count.nf' params(hts_count_type:params.hts_count_type,
                                                                    hts_group_features:params.hts_group_features,
                                                                    optional:params.count.toolOptions, 
                                                                    singleEnd:params.singleEnd, 
                                                                    stranded:params.stranded, 
                                                                    unstranded:params.unstranded, 
                                                                    revstranded:params.revstranded)
include AlignReads from './NextflowModules/STAR/2.7.3a/AlignReads.nf' params(singleEnd:params.singleEnd, 
									                                                           optional:params.star.toolOptions)
include Index from './NextflowModules/Sambamba/0.6.8/Index.nf' params(params)
include QuantMerge from './NextflowModules/Salmon/1.2.1/QuantMerge.nf' params( optional: params.salmon_quantmerge.toolOptions )
include Quant from './NextflowModules/Salmon/1.2.1/Quant.nf' params(singleEnd: params.singleEnd,
                                                                     stranded: params.stranded,
                                                                     unstranded: params.unstranded,
                                                                     revstranded: params.revstranded,
                                                                     saveUnaligned: params.saveUnaligned,
								                                                     optional: params.salmon_quant.toolOptions )
                                                                  

include MergeFastqLanes from './NextflowModules/Utils/MergeFastqLanes.nf' params(params)
include mergeHtseqCounts from './utils/mergeHtseqCounts.nf' params(params)
include EdgerNormalize as fc_norm from './utils/bioconductor/edger/3.28.0/normalize.nf' params( tool:"fc" )
include FeatureCounts from './NextflowModules/Subread/2.0.0/FeatureCounts.nf' params( optional:params.fc.toolOptions,
										                                                                  biotypeQC:params.biotypeQC,
                                                                                      singleEnd: params.singleEnd,
                                                                                      stranded: params.stranded,
                                                                                      unstranded: params.unstranded,
                                                                                      revstranded: params.revstranded,
                                                                                      fc_group_features: params.fc_group_features,
                                                                                      fc_count_type: params.fc_count_type,
                                                                                      fc_group_features_type: params.fc_group_features_type,
                                                                                      fc_extra_attributes : params.fc_extra_attributes, 
                                                                                      gencode: params.gencode)
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

workflow {
  main :  
    run_name = params.fastq_path.split('/')[-1]
    fastq_files = extractAllFastqFromDir(params.fastq_path).map { [it[0],it[1],it[4]]}
    //Get necessary files
    genome_gtf = Channel
        .fromPath(params.genome_gtf, checkIfExists: true)
        .ifEmpty { exit 1, "GTF file not found: ${params.genome_gtf}"}
    genome_fasta = Channel
        .fromPath(params.genome_fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Fasta file not found: ${params.genome_fasta}"}
    //Determine required inputs
    if (params.star_index && params.runMapping) {
      star_index = Channel
            .fromPath(params.star_index, checkIfExists: true)
            .ifEmpty { exit 1, "STAR index not found: ${params.star_index}"}
    } else if (!params.star_index && params.runMapping) {
      //Create STAR Index
      GenomeGenerate ( genome_fasta, genome_gtf )
      star_index = GenomeGenerate.out.star_index
    }
    if (params.genome_bed && params.runPostQC && params.runMapping) {
      //Create bed12 index file
      genome_bed = Channel
            .fromPath(params.genome_bed, checkIfExists: true)
            .ifEmpty { exit 1, "Bed12 file not found: ${params.genome_bed}"}
    } else if ( !params.genome_bed && params.runPostQC && params.runMapping) {
        GtfToGenePred ( genome_gtf)
        GenePredToBed ( GtfToGenePred.out.genome_genepred )
        genome_bed = GenePredToBed.out.genome_bed12
    }
    if ( params.salmon_index && params.runSalmon) {
       salmon_index = Channel
            .fromPath(params.salmon_index, checkIfExists: true)
            .ifEmpty { exit 1, "Transcripts fasta not found: ${params.salmon_index}"}
    } else if ( !params.salmon_index && params.runSalmon ) {
        transcripts_fasta = Channel
            .fromPath(params.transcripts_fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Fasta file not found: ${params.transcripts_fasta}"}
        SalmonIndex ( transcripts_fasta )
        salmon_index = SalmonIndex.out.salmon_index
    }
    if (params.scatter_interval_list && params.runGATK4_HC ) {
      scatter_interval_list = Channel
        .fromPath( params.scatter_interval_list, checkIfExists: true)
        .ifEmpty { exit 1, "Scatter intervals not found: ${params.scatter_interval_list}"}
    } else if ( !params.scatter_interval_list && params.runGATK4_HC ) {
        genome_dict = Channel
              .fromPath( params.genome_dict, checkIfExists: true)
              .ifEmpty { exit 1, "Genome dictionary not found: ${params.genome_dict}"}
        genome_index = Channel
              .fromPath(params.genome_fasta + '.fai', checkIfExists: true)
              .ifEmpty { exit 1, "Fai file not found: ${params.genome_fasta}.fai"}
        CreateIntervalList( genome_index, genome_dict )
        scatter_interval_list = CreateIntervalList.out.genome_interval_list
    }
    if ( params.runSortMeRna) {
        rRNA_database = file(params.rRNA_database_manifest)
        //if (rRNA_database.isEmpty()) {exit 1, "File ${rRNA_database.getName()} is empty!"}
        sortmerna_fasta = Channel
            .from( rRNA_database.readLines() )
            .map { row -> file(row) }
    }
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
    // Determine final fastqs files
    if ( params.runTrimGalore && params.runSortMeRna ) {
      TrimGalore(fastq_files) 
      SortMeRNA(TrimGalore.out.trimmed_fastqs, sortmerna_fasta.collect())
      final_fastqs = SortMeRNA.out.non_rRNA_fastqs

    } else if ( params.runTrimGalore && !params.runSortMeRna ) {
        TrimGalore(fastq_files)
        final_fastqs = TrimGalore.out.trimmed_fastqs 

    } else if ( !params.runTrimGalore &&  params.runSortMeRna ) {
        SortMeRNA(fastq_files, sortmerna_fasta.collect() )
        final_fastqs = SortMeRNA.out.non_rRNA_fastqs 

    } else {
        final_fastqs = fastq_files
    }
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
    if ( params.runMapping ) {
      AlignReads( fastqs_transformed, star_index.collect(), genome_gtf.collect() )
      Index(AlignReads.out.star_aligned.map {sample_id, rg_id, bam ->
                                             [sample_id, bam] })
      mapped = AlignReads.out.star_aligned.join(Index.out)
    }
    if ( params.runPostQC) {
      if (params.runMapping) {
        post_mapping_QC(mapped.map { sample_id, rg_id, bam, bai -> 
                                    [sample_id, bam, bai] }, 
                                    genome_bed.collect())
      } else {
          exit 1, "PostQC requires alignment step. Please enable runMapping!"
      } 
    } 
    if ( params.runHTSeqCount) {
      if ( params.runMapping) {
        Count(mapped.map { sample_id, rg_id, bam, bai -> 
                          [sample_id, bam, bai] }, 
                          genome_gtf.collect())
        mergeHtseqCounts( run_name, Count.out.hts_counts_raw.collect() )
      } else {
          exit 1, "htseq-count requires alignment step. Please enable runMapping!"
      }
    } 
    if ( params.runFeatureCounts) {
      if (params.runMapping) {
        FeatureCounts(run_name, mapped.map { it[2] }.collect(), genome_gtf.collect())
        if ( params.normalize_counts ) {
          fc_norm( run_name, FeatureCounts.out.fc_raw )
        }
      } else {
          exit 1, "featureCounts requires alignment step. Please enable runMapping!"
      } 
    } 
    if ( params.runMarkDup) {
      if (params.runMapping) {
        markdup_mapping(mapped)
      } else {
          exit 1, "Markdup requires alignment step. Please enable runMapping!"
      } 
    }         
    if ( params.runSalmon ) {
      Quant ( MergeFastqLanes (fastqs_transformed ), salmon_index.collect() )
      QuantMerge ( run_name, Quant.out.map { it[1] }.collect() )
    }
    if ( params.runGATK4_HC ) {
      if (params.runMapping && params.runMarkDup) {
        SplitIntervals( 'no-break', scatter_interval_list)
        SplitNCigarReads(markdup_mapping.out)
        if ( params.runGATK4_BQSR) {
          //Perform BSQR
          gatk4_bqsr(SplitNCigarReads.out, SplitIntervals.out.flatten())
          gatk4_hc(gatk4_bqsr.out[0], SplitIntervals.out.flatten(), run_name)
        } else {
            gatk4_hc(SplitNCigarReads.out, SplitIntervals.out.flatten(), run_name)
        }      
      }  else {
       exit 1, "GATK4 requires alignment, markdup step. Please enable runMapping and runMarkDup!"
      }     
    }  
    if ( params.runMultiQC ) {
      //Create empty Channels for optional steps
      trim_logs = Channel.empty()
      fastqc_logs = Channel.empty()
      sortmerna_logs = Channel.empty()
      star_logs = Channel.empty()
      post_qc_logs = Channel.empty()
      hts_logs = Channel.empty()
      fc_logs = Channel.empty()
      salmon_logs = Channel.empty()
      //Get options
      if (  params.runTrimGalore) {
        trim_logs = TrimGalore.out.trimming_report
        fastqc_logs = TrimGalore.out.trimming_fastqc
      }
      if (  params.runSortMeRna ) {
        //Currently not working with MultiQc 1.8
        //sortmerna_logs = SortMeRNA.out.sortmerna_report
      }
      if (  params.runMapping ) {
        star_logs =  AlignReads.out.star_log.mix(AlignReads.out.star_final_log)
      }
      if ( params.runHTSeqCount &&  params.runMapping) {
        hts_logs = Count.out.hts_counts_raw
      }
      if (  params.runFeatureCounts &&  params.runMapping ) {
        fc_logs = FeatureCounts.out.fc_summary
      }
      if ( params.runPostQC && params.runMapping ) {
        post_qc_logs = post_mapping_QC.out[1].map { it[1] }.mix(post_mapping_QC.out[0].map { it[1] })
      }
      if ( params.runSalmon) {
        salmon_logs = Quant.out.salmon_quants
      }
      multiqc_report( run_name,
		      fastqc_logs,
                      trim_logs,
                      sortmerna_logs,
                      star_logs,
                      post_qc_logs,
                      hts_logs,
                      fc_logs,
                      salmon_logs ) 		      
     }

}
