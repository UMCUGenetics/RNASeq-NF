#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*===========================================================
                UMCUGenetics + EMC / RNASeq-NF
===========================================================
#### Homepage / Documentation
https://github.com/UMCUGenetics/RNASeq-NF
----------------------------------------------------------------------------------------
*/
def helpMessage() {
    // Log colors ANSI codes
    c_reset = "\033[0m";
    c_dim = "\033[2m";
    c_black = "\033[0;30m";
    c_green = "\033[0;32m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_purple = "\033[0;35m";
    c_white = "\033[0;37m";

    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run UMCUGenetics/RNASeq-NF --fastq_path <fastq_dir> --out_dir <output_dir> --genome_config <path/to/genome.config
${c_blue}    Mandatory arguments: ${c_reset}
${c_yellow}        --fastq_path [str] ${c_reset}              Path to a directory containing fastq files.
                                                    Files should be named in the following format: xxx_xxx_xxxx
${c_yellow}        --out_dir [str] ${c_reset}                 The output directory where the results will be saved
${c_blue}    Standard options: ${c_reset}
         --profile [str]                 Configuration profile to use, leave empty to run locally.
                                              Available: slurm, SGE, singularity.
         --genome_config [path]          Path to genome configuration file containing options from ${c_blue}standard references${c_reset}.
         --singleEnd [bool]             Specifies that the input is from single-end experiment(s). (Default: false)
         --unstranded [bool]             Specifies that the input is from an unstranded library prep. (Default: true)
         --stranded [bool]               Specifies that the input is from an forward-stranded library prep. (Default: false)
         --revstranded [bool]            Specifies that the input is from an reverse-stranded library prep. (Default: false)
${c_blue}    Standard references: ${c_reset}
      If not specified in the configuration file or you wish to overwrite any of standard references.
${c_yellow}        --genome_fasta [path] ${c_reset}           Path to genome sequence file (FASTA).
${c_yellow}        --genome_gtf [path] ${c_reset}             Path to GTF file containing genomic annotations.
${c_yellow}        --genome_bed [path] ${c_reset}             Path to BED12-format of the supplied GTF (auto-generated from supplied GTF if not given).
        --star_index [path]              Path to STAR index (generated automatically if not given).
        --gencode [bool]                 Specifies if the supplied GTF is from GENCODE. (Default: false)
${c_blue}    FastQC: ${c_reset}
      Perform FastQC on the unaligned sequencing reads before and, optionally, after trimming.
${c_green}        --runFastQC [bool] ${c_reset}          Run FastQC. (Default: true)
        --options.FastQC [str]       Additional custom options given to FastQC.
${c_blue}    TrimGalore: ${c_reset}
      Trims sequence adapters from sequencing reads and filters low-complexity and small reads.
${c_green}        --runTrimGalore [bool] ${c_reset}          Run TrimGalore. (Default: true)
        --options.TrimGalore [str]       Additional custom options given to TrimGalore.
${c_blue}    SortMeRNA: ${c_reset}
      Removes sequencing reads from ribosomal origins.
${c_green}        --runSortMeRNA [bool] ${c_reset}           Run SortMeRNA. (Default: true)
${c_yellow}        --rRNA_database_manifest [path]${c_reset}  Path to rRNA database files.
        --options.SortMeRNA [str]        Additional custom options given to SortMeRNA.
${c_blue}    Alignment - STAR/MarkDup: ${c_reset}
      Performs alignment of sequencing reads against the genome using STAR.
${c_green}        --runMapping [bool] ${c_reset}                Run STAR. (Default: true)
        --options.STAR [str]             Additional custom options given to STAR.
${c_blue}    Post-alignment QC: ${c_reset}
      Various QC to perform after alignment to assess quality of alignment and sequencing read input.
${c_green}        --runRSeQC_TIN [bool] ${c_reset}                 Run tin.py to assess distribution of reads over gene-body. (Default: true)
${c_green}        --runPostQC [bool] ${c_reset}                  Run RSeQC components: inner_distance, read_distribution, infer_experiment, junction_annotation, bam_stat, junction_saturation and read_duplication. Run preseq to predict and estimate the complexity of the sequencing library (Default: true)
${c_blue}    Counting - SubRead / FeatureCounts: ${c_reset}
      Read counting, per BAM file, per <fc_group_features> by counting all reads overlapping with <fc_count_type>.
${c_green}        --runFeatureCounts [bool] ${c_reset}       Run FeatureCounts. (Default: true)
        --fc_group_features [str]        Feature to summarize reads on. (Default: gene_id)
        --fc_count_type [str]            Feature to count overlapping reads, and subsequently summarized by --fc_group_features. (Default: exon)
        --fc_group_features_type [str]   GTF biotype field for subread featureCounts (Default: gene_biotype)
        --normalize_counts [bool]        Enable edgeR RPKM/CPM normalization for featureCounts (Default: true) 
        --options.FeatureCounts [str]    Additional custom options given to FeatureCounts.
${c_blue}    Salmon: ${c_reset}
      Performs transcript alignment and quantification of the expression of transcripts, per isoform.
${c_green}        --runSalmon [bool] ${c_reset}              Run Salmon. (Default: false)
${c_yellow}        --transcripts_fasta [path] ${c_reset}      Path to transcripts in FASTA format.
        --salmon_index [path]            Path to Salmon Index (auto-generated if not given).
      Additional custom options given to Salmon submodules.
        --options.Salmon_quant [str]
        --options.Salmon_index [str] 
        --options.Salmon_quantmerge [str]
${c_blue}    GATK (v4) - Germline variant calling: ${c_reset}
      Performs germline variant calling using the RNA-Seq best-practices as established by GATK.
${c_green}        --runGermlineCallingGATK [bool] ${c_reset} Run GATK4 for (germline) variant calling. (Default: false)
${c_yellow}        --scatter_interval_list [path] ${c_reset}  Path to scatter.interval_list (required for GATK4)
        --genome_known_sites [path]      Path to snp_sites.vcf (optional for use in GATK4 BQSR)
      Additional custom options given to GATK4 tools.
          --options.GATK4_SplitIntervals [str]
          --options.GATK4_HaplotypeCaller [str]
          --options.GATK4_VariantFiltration [str]
${c_blue}    GATK (v4) - Base quality score recalibration (BQSR): ${c_reset}
      Performs BQSR.
${c_green}        --runGATK4_BQRS [bool] ${c_reset}  Run BQRS to recalibrate base quality scores. (Default: false)
        --options.GATK4_BQRS [str]              Additional custom options given to BQRS.
${c_blue}    MultiQC: ${c_reset}
      Generate a MultiQC report which combined various QC reports into a single report.
${c_green}        --runMultiQC [bool] ${c_reset}             Perform MultiQC to generate a single report containing various QC logs.
        --options.MultiQC [str]          Additional custom options given to MultiQC.
    """.stripIndent()
}


/*=================================
          Input validation
=================================*/

// Show help message and exit.
if(params.help){
  helpMessage()
  exit 0
}

// Minimal required parameters.
if (!params.out_dir) {
   exit 1, "Output directory not found, please provide the correct path! (--out_dir)"
}

if (!params.fastq_path) {
  exit 1, "fastq files not found, please provide the correct path! (--fastq_path)"
}

if (!params.genome_fasta) {
  exit 1, "Genome fasta not found, please provide the correct path! (--genome_fasta)"
} else {
  // Try importing.
  genome_fasta = Channel
      .fromPath(params.genome_fasta, checkIfExists: true)
      .ifEmpty { exit 1, "Fasta file not found: ${params.genome_fasta}"}
}
if (!params.genome_gtf) {
  exit 1, "Genome GTF not found, please provide the correct path! (--genome_gtf)"
} else {
  // Try importing.
  genome_gtf = Channel
      .fromPath(params.genome_gtf, checkIfExists: true)
      .ifEmpty { exit 1, "GTF file not found: ${params.genome_gtf}"}
}

//Start workflow
workflow {
  main :
    //Set run and retrieve input fastqs
    include extractAllFastqFromDir from './NextflowModules/Utils/fastq.nf' params(params)  
    run_name = params.fastq_path.split('/')[-1]
    fastq_files = extractAllFastqFromDir(params.fastq_path).map { [it[0],it[1],it[4]]}
   //Pipeline log info
    params.version = "1.0.0-rc1"
    log.info """=======================================================
    RNASeq-NF ${params.version}"
    ======================================================="""
    def summary = [:]
    summary['Pipeline Name']  = 'RNASeq-NF'
    summary['Pipeline Version'] = params.version
    summary['Run Name']     = run_name
    summary['Fastq dir']   = params.fastq_path
    summary['Genome fasta']   = params.genome_fasta
    summary['Genome GTF']   = params.genome_gtf
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
    include pre_processing from './sub-workflows/pre_processing.nf' params(params)
    pre_processing ( fastq_files )
    //Logs
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
        include markdup_mapping from './sub-workflows/mapping_deduplication.nf' params(params)
        mapped = markdup_mapping( fastqs_transformed, genome_fasta, genome_gtf )
        star_logs = mapped.logs
    } else {
        star_logs = Channel.empty()
    }
    // # 3) Post-mapping QC
    if ( params.runPostQC) {
      if (params.runMapping) {
          include post_mapping_QC from './sub-workflows/post_mapping_QC.nf' params(params)
          post_mapping_QC( mapped.bam_sorted.map { sample_id, rg_id, bam, bai -> [sample_id, bam, bai] }, genome_gtf )
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
         include alignment_based_quant from './sub-workflows/alignment_based_quant.nf' params(params)
         alignment_based_quant ( run_name, mapped.bam_sorted.map { it[2] }, genome_gtf )
         fc_logs = alignment_based_quant.out.fc_summary
        } else {
          exit 1, "featureCounts requires alignment step. Please enable runMapping!"
      } 
    } else {
          fc_logs = Channel.empty()
    }
    // # 5) Salmon    
    if ( params.runSalmon ) {
      include alignment_free_quant from './sub-workflows/alignment_free_quant.nf' params(params)
      alignment_free_quant( fastqs_transformed, run_name )
      salmon_logs = alignment_free_quant.out.logs
    } else {
        salmon_logs = Channel.empty()
    }
    // # 6) GATK4 germline variant calling with optional BQSR
    if ( params.runGermlineCallingGATK ) {
      if ( params.runMapping ) {
          include gatk_germline_calling from './sub-workflows/gatk_germline_calling.nf' params(params)
          gatk_germline_calling( run_name, mapped.bam_dedup.map {sample_id, rg_id, bam, bai -> 
                                                                 [sample_id, bam, bai] }) 
      } else {
            exit 1, "GATK4 requires alignment, markdup step. Please enable runMapping and runMarkDup!"
      }        
    }
    // # 7) MultiQC report  
    if ( params.runMultiQC ) {
        include multiqc_report from './sub-workflows/multiqc_report.nf' params(params)
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
