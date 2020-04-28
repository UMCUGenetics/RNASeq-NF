#!/usr/bin/env nextflow

nextflow.preview.dsl=2
include TrimGalore from './NextflowModules/TrimGalore/0.6.5/TrimGalore.nf' params( optional: params.trimgalore.toolOptions, 
                                                                                   singleEnd: params.singleEnd )
include GenomeGenerate from './NextflowModules/STAR/2.7.3a/GenomeGenerate.nf' params(params)
include Index as SalmonIndex from './NextflowModules/Salmon/0.13.1/Index.nf' params( gencode: params.gencode )
include GtfToGenePred from './NextflowModules/ucsc/377/gtfToGenePred/GtfToGenePred.nf' params(params)
include GenePredToBed from './NextflowModules/ucsc/377/genePredToBed/GenePredToBed.nf' params(params)
include CreateSequenceDictionary from './NextflowModules/Picard/2.22.0/CreateSequenceDictionary.nf' params(params)
include CreateIntervalList from './NextflowModules/Utils/CreateIntervaList.nf' params(params)
include getExonLenghts from './utils/getExonLengths.nf' params(params)
include extractAllFastqFromDir from './NextflowModules/Utils/fastq.nf' params(params)
include post_mapping_QC from './sub-workflows/post_mapping_QC.nf' params(params)
include markdup_mapping from './sub-workflows/mapping_deduplication.nf' params(params)
include multiqc_report from './sub-workflows/multiqc_report.nf' params(params)
include SortMeRna from './NextflowModules/SortMeRNA/4.2.0/SortMeRna.nf' params(singleEnd:params.singleEnd)
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
include Quant from './NextflowModules/Salmon/0.13.1/Quant.nf' params(singleEnd: params.singleEnd,
                                                                     stranded: params.stranded,
                                                                     unstranded: params.unstranded,
                                                                     revstranded: params.revstranded,
                                                                     saveUnaligned: params.saveUnaligned)

include mergeFastqLanes from './NextflowModules/Utils/mergeFastqLanes.nf' params(params)
include mergeHtseqCounts from './utils/mergeHtseqCounts.nf' params(params)
include mergeSalmonCounts from './utils/mergeSalmonCounts.nf' params(params)
include rpkm as hts_rpkm from './utils/bioconductor/edger/3.28.0/rpkm.nf' params(tool:"hts")
include rpkm as fc_rpkm from './utils/bioconductor/edger/3.28.0/rpkm.nf' params(tool:"fc")
include FeatureCounts from './NextflowModules/subread/2.0.0/FeatureCounts.nf' params( optional:params.fc.toolOptions,
                                                                                      singleEnd: params.singleEnd,
                                                                                      stranded: params.stranded,
                                                                                      unstranded: params.unstranded,
                                                                                      revstranded: params.revstranded,
                                                                                      fc_group_features: params.fc_group_features,
                                                                                      fc_count_type: params.fc_count_type,
                                                                                      fc_group_features_type: params.fc_group_features_type,
                                                                                      fc_extra_attributes : params.fc_extra_attributes, 
                                                                                      gencode: params.gencode)



if (!params.out_dir) {
   exit 1, "Output directory not found. Please provide the correct path!"
}

if (!params.fastq_path) {
  exit 1, "fastq files not found. Please provide the correct path!"
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
    genome_index = Channel
        .fromPath(params.genome_fasta + '.fai', checkIfExists: true)
        .ifEmpty { exit 1, "Fai file not found: ${params.genome_fasta}.fai"}

    if (params.gene_len && params.norm_rpkm) {
      exon_lengths = params.gene_len
    } else if (!params.gene_len && params.norm_rpkm ) {
      getExonLenghts( genome_gtf)
      exon_lengths = getExonLenghts.out
    } 
    if (params.star_index && !params.skipMapping) {
      star_index = Channel
            .fromPath(params.star_index, checkIfExists: true)
            .ifEmpty { exit 1, "STAR index not found: ${params.star_index}"}
    } else if (!params.star_index && !params.skipMapping) {
      //Create STAR Index
      GenomeGenerate ( genome_fasta, genome_gtf )
      star_index = GenomeGenerate.out
    }
    if (params.genome_bed && !params.skipPostQC && !params.skipMapping) {
      //Create bed12 index file
      genome_bed = Channel
            .fromPath(params.genome_bed, checkIfExists: true)
            .ifEmpty { exit 1, "Bed12 file not found: ${params.genome_bed}"}
    } else if ( !params.genome_bed && !params.skipPostQC && !params.skipMapping) {
        GtfToGenePred ( genome_gtf)
        GenePredToBed ( GtfToGenePred.out )
        genome_bed = GenePredToBed.out
    }
    if ( params.salmon_index && !params.skipSalmon) {
       salmon_index = Channel
            .fromPath(params.salmon_index, checkIfExists: true)
            .ifEmpty { exit 1, "Transcripts fasta not found: ${params.salmon_index}"}
    } else if ( !params.salmon_index && !params.skipSalmon ) {
        transcripts_fasta = Channel
            .fromPath(params.transcripts_fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Fasta file not found: ${params.transcripts_fasta}"}
      SalmonIndex ( transcripts_fasta )
      salmon_index = SalmonIndex.out
    }
    if (params.scatter_interval_list && !params.skipGATK4_HC ) {
      scatter_interval_list = Channel
        .fromPath( params.scatter_interval_list, checkIfExists: true)
        .ifEmpty { exit 1, "Scatter intervals not found: ${params.scatter_interval_list}"}
    } else if ( !params.scatter_interval_list && !params.skipGATK4_HC) {
        CreateSequenceDictionary (genome_fasta)
        CreateIntervalList(genome_index, CreateSequenceDictionary.out )
        scatter_interval_list = CreateIntervalList.out
    }
    if ( !params.skipSortMeRna) {
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
    //Log channels
    trim_logs = Channel.empty()
    fastqc_logs = Channel.empty()
    sortmerna_logs = Channel.empty()
    // Determine final fastqs files
    if ( !params.skipTrimGalore && !params.skipSortMeRna ) {
      TrimGalore(fastq_files) 
      SortMeRna(TrimGalore.out.map{ sample_id, rg_id, reads, log, fqc_report -> [sample_id, rg_id, reads] }, 
                sortmerna_fasta.collect() )
      final_fastqs = SortMeRna.out.map{ [it[0],it[1],it[2]] }

    } else if ( !params.skipTrimGalore && params.skipSortMeRna ) {
        TrimGalore(fastq_files)
        final_fastqs = TrimGalore.out.map{ sample_id, rg_id, reads, log, fqc_report -> [sample_id, rg_id, reads] }

    } else if ( params.skipTrimGalore &&  !params.skipSortMeRna ) {
        SortMeRna(fastq_files, 
                  sortmerna_fasta.collect() )
        final_fastqs = SortMeRna.out.map{ [it[0],it[1],it[2]] }

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
    if (!params.skipMapping) {
      AlignReads( fastqs_transformed, star_index.collect(), genome_gtf.collect() )
      Index(AlignReads.out.map { sample_id, bams, unmapped, log1, log2, tab -> [sample_id, bams] })
      mapped = AlignReads.out.join(Index.out)
    }
    if (!params.skipPostQC && !params.skipMapping) {
      post_mapping_QC(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, bams, bai] }, genome_bed.collect())
    }
    if (!params.skipHTSeqCount && !params.skipMapping) {
      Count(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, bams, bai] }, genome_gtf.collect())
      mergeHtseqCounts( run_name, Count.out.map { it[1] }.collect())
      if ( params.norm_rpkm) {
        hts_rpkm( run_name, mergeHtseqCounts.out, exon_lengths)
      }
    }
    if (!params.skipFeatureCounts && !params.skipMapping) {
      FeatureCounts(run_name, AlignReads.out.map { it[1] }.collect(), genome_gtf.collect())
      //if ( params.norm_rpkm ) {
      //  fc_rpkm( run_name, FeatureCounts.out.map { it[1] }, exon_lengths)
      //}
    }
    if (!params.skipMarkDup && !params.skipMapping) {
      markdup_mapping(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, sample_id, bams, bai] })
    }
    if (!params.skipSalmon) {
      Quant ( mergeFastqLanes (fastqs_transformed ), salmon_index.collect() ) 
      mergeSalmonCounts ( run_name, Quant.out.map { it[1] }.collect())
    }
    if (!params.skipMapping && !params.skipMarkDup && !params.skipGATK4_HC) {
          SplitIntervals( 'no-break', scatter_interval_list)
          SplitNCigarReads(markdup_mapping.out)
          if (!params.skipGATK4_BQSR) {
            //Perform BSQR
            gatk4_bqsr(SplitNCigarReads.out, SplitIntervals.out.flatten())
            gatk4_hc(gatk4_bqsr.out[0], SplitIntervals.out.flatten())
          } else {
              gatk4_hc(SplitNCigarReads.out, SplitIntervals.out.flatten())
          }      
    }
    if ( !params.skipMultiQC ) {
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
      if ( !params.skipTrimGalore) {
        trim_logs = TrimGalore.out.map { it[3] }
        fastqc_logs = TrimGalore.out.map { it[4] }
      }
      if ( !params.skipSortMeRna ) {
        //Currently not working with MultiQc 1.8
        //sortmerna_logs = SortMeRna.out.map { it[3] }
      }
      if ( !params.skipMapping) {
        star_logs =  AlignReads.out.map{ [it[3], it[4]] }
      }
      if ( !params.skipHTSeqCount && !params.skipMapping) {
        hts_logs = Count.out.map { it[1] }
      }
      if ( !params.skipFeatureCounts && !params.skipMapping) {
        fc_logs = FeatureCounts.out.map { it[4]}
      }
      if ( !params.skipPostQC && !params.skipMapping ) {
        post_qc_logs =  post_mapping_QC.out[1].map { it[1] }.mix(post_mapping_QC.out[0].map { it[1] })
      }
      if ( !params.skipSalmon) {
        salmon_logs = Quant.out.map { it[1] }
      }
      multiqc_report( fastqc_logs,
                      trim_logs,
                      sortmerna_logs,
                      star_logs,
                      post_qc_logs,
                      hts_logs,
                      fc_logs,
                      salmon_logs ) 		      
    }

}
