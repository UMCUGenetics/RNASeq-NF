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
         --single_end [bool]             Specifies that the input is from single-end experiment(s). (Default: false)
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
${c_green}        --runFastQCTrimmed [bool] ${c_reset}   Run FastQC after trimming (if enabled). (Default: true)
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

${c_blue}    Alignment - STAR: ${c_reset}
      Performs alignment of sequencing reads against the genome using STAR.
${c_green}        --runSTAR [bool] ${c_reset}                Run STAR. (Default: true)
${c_green}        --runMarkDup [bool] ${c_reset}             Runs Sambamba MarkDuplicates. (Default: true)
        --options.STAR [str]             Additional custom options given to STAR.

${c_blue}    Post-alignment QC: ${c_reset}
      Various QC to perform after alignment to assess quality of alignment and sequencing read input.
${c_green}        --runTIN [bool] ${c_reset}                 Run tin.py to assess distribution of reads over gene-body. (Default: true)
${c_green}        --RSeQC [bool] ${c_reset}                  Run RSeQC components: inner_distance, read_distribution, infer_experiment, junction_annotation, bam_stat, junction_saturation and read_duplication. (Default: true)
${c_green}        --runPreseq [bool] ${c_reset}              Run preseq to predict and estimate the complexity of the sequencing library. (Default: true)

${c_blue}    Counting - SubRead / FeatureCounts: ${c_reset}
      Read counting, per BAM file, per <fc_group_features> by counting all reads overlapping with <fc_count_type>.
${c_green}        --runFeatureCounts [bool] ${c_reset}       Run FeatureCounts. (Default: true)
        --fc_group_features [str]        Feature to summarize reads on. (Default: gene_id)
        --fc_count_type [str]            Feature to count overlapping reads, and subsequently summarized by --fc_group_features. (Default: exon)
        --fc_group_features_type [str]   GTF biotype field for subread featureCounts (Default: gene_biotype)
        --normalize_counts [bool]        Enable edgeR RPKM/CPM normalization for featureCounts (Default: true) ????
        --options.FeatureCounts [str]    Additional custom options given to FeatureCounts.

${c_blue}    Salmon: ${c_reset}
      Performs transcript alignment and quantification of the expression of transcripts, per isoform.
${c_green}        --runSalmon [bool] ${c_reset}              Run Salmon. (Default: false)
${c_yellow}        --transcripts_fasta [path] ${c_reset}      Path to transcripts in FASTA format.
        --salmon_index [path]            Path to Salmon Index (auto-generated if not given).
        --options.Salmon [str]           Additional custom options given to Salmon.

${c_blue}    GATK (v4) - Germline variant calling: ${c_reset}
      Performs germline variant calling using the RNA-Seq best-practices as established by GATK.
${c_green}        --runGermlineCallingGATK [bool] ${c_reset}  Run GATK4 for (germline) variant calling. (Default: false)
${c_yellow}        --scatter_interval_list [path] ${c_reset}  Path to scatter.interval_list (required for GATK4)
        --genome_known_sites [path]      Path to snp_sites.vcf (optional for use in GATK4 BQSR)
        --options.GATKGermline [str]             Additional custom options given to GATK4.

${c_blue}    GATK (v4) - Base quality score recalibration (BQSR): ${c_reset}
      Performs BQSR.
${c_green}        --runBQRS [bool] ${c_reset}  Run BQRS to recalibrate base quality scores. (Default: false)

${c_blue}    MultiQC: ${c_reset}
      Generate a MultiQC report which combined various QC reports into a single report.
${c_green}        --runMultiQC [bool] ${c_reset}             Perform MultiQC to generate a single report containing various QC logs.

${c_blue}    Other options: ${c_reset}
      --email [email]                     Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits.
      --email_on_fail [email]             Same as --email, except only send mail if the workflow is not successful.
      --max_multiqc_email_size [str]      Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      --name [str]                        Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
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

// Minimally required parameters.
if (!params.out_dir) {
   exit 1, "Output directory not found, please provide the correct path! (--out_dir)"
}

if (!params.fastq_path) {
  exit 1, "fastq files not found, please provide the correct path! (--fastq_path)"
}

if (!params.genome_fasta) {
  exit 1, "Genome fasta not found, please provide the correct path! (--genome_fasta)"
}else{

  // Try importing.
  genome_fasta = Channel
      .fromPath(params.genome_fasta, checkIfExists: true)
      .ifEmpty { exit 1, "Fasta file not found: ${params.genome_fasta}"}
}

if (!params.genome_gtf) {
  exit 1, "Genome GTF not found, please provide the correct path! (--genome_gtf)"
}else{

  // Try importing.
  genome_gtf = Channel
      .fromPath(params.genome_gtf, checkIfExists: true)
      .ifEmpty { exit 1, "GTF file not found: ${params.genome_gtf}"}
}

// Run workflow.
workflow {
  main :

    // Include required modules.
    include extractAllFastqFromDir from './NextflowModules/Utils/fastq.nf' params(params)

    // Determine run name and retrieve all fastq files from input directory.
    run_name = params.fastq_path.split('/')[-1]
    fastq_files = extractAllFastqFromDir(params.fastq_path).map { [it[0],it[1],it[4]]}

    //Determine required inputs.




    // Output general information of run.
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

    final_fastqs = fastq_files

    // ToDo - Make separate NF module.
    // FastQC - Summary of raw reads.
    if ( params.runFastQC ) {
      // Include required modules.


      // Run workflow.


      // Logging
      // logs_FastQC = FastQC.out.map { it[4] }

    }else{

      // Create empty log.
      logs_FastQC = Channel.empty()
    }


    // Trimming reads - Trimmomatic
    if ( params.runTrimGalore ) {

      // Include required modules.
      include TrimGalore from './NextflowModules/TrimGalore/0.6.5/TrimGalore.nf' params(
        optional: params.trimgalore.toolOptions,
        singleEnd: params.singleEnd
      )

      // Run workflow.
      TrimGalore(fastq_files)
      final_fastqs = TrimGalore.out.map{ sample_id, rg_id, reads, log, fqc_report -> [sample_id, rg_id, reads] }

      // Logging
      logs_TrimGalore = TrimGalore.out.map { it[3] }

    }else{

      // Create empty log.
      logs_TrimGalore = Channel.empty()
    }


    // Removal of ribosomal RNA - SortMeRNA.
    if(params.runSortMeRNA){

      if (params.rRNA_database) {

        rRNA_database = file(params.rRNA_database_manifest)
        if (rRNA_database.isEmpty()) {exit 1, "File ${rRNA_database.getName()} is empty!"}

        sortmerna_fasta = Channel
            .from( rRNA_database.readLines() )
            .map { row -> file(row) }

      }else{
        exit 1, "Path to Ribosomal RNA database(s) not given, please provide the correct parameter! (--rRNA_database)"
      }

      // Include required modules.
      include SortMeRna from './NextflowModules/SortMeRNA/4.2.0/SortMeRna.nf' params(
        singleEnd:params.singleEnd
      )

      // Run workflow.
      SortMeRna(fastq_files, sortmerna_fasta.collect() )
      final_fastqs = SortMeRNA.out.map{ [it[0],it[1],it[2]] }

      // Logging.
      logs_SortMeRNA = SortMeRna.out.map { it[3] }

    }else{

      // Create empty log.
      logs_SortMeRNA = Channel.empty()

    }


    // Transform output channels
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


    // Alignment using STAR.
    if ( params.runSTAR ) {

      // Include required modules.
      include GenomeGenerate as GenomeGenerateSTAR from './NextflowModules/STAR/2.7.3a/GenomeGenerate.nf' params(params)
      include AlignReads as AlignReadsSTAR from './NextflowModules/STAR/2.7.3a/AlignReads.nf' params(
        singleEnd:params.singleEnd,
        optional:params.star.toolOptions
      )

      // Check / Generate STAR Index
      if (params.star_index) {

        star_index = Channel
              .fromPath(params.star_index, checkIfExists: true)
              .ifEmpty { exit 1, "Supplied STAR index not found: ${params.star_index} (--star_index)"}

      } else if (!params.star_index) {

        include Index as IndexSambamba from './NextflowModules/Sambamba/0.7.1/Index.nf' params(params)

        //Create STAR Index
        GenomeGenerateSTAR ( genome_fasta, genome_gtf )
        star_index = GenomeGenerateSTAR.out

      }

      // Run workflow.
      AlignReadsSTAR( fastqs_transformed, star_index.collect(), genome_gtf.collect() )
      IndexSambamba(AlignReadsSTAR.out.map { sample_id, bams, unmapped, log1, log2, tab -> [sample_id, bams] })
      mapped = AlignReadsSTAR.out.join(IndexSambamba.out)

      // Logging.
      logs_STAR = AlignReadsSTAR.out.map{ [it[3], it[4]] }

    }else{

      // Create empty log.
      logs_STAR = Channel.empty()

    }


    // Marking duplicates - Sambamba.
    if ( params.runMarkDup) {

      // Input validation.
      if (!params.runSTAR){
        exit 1, "MarkDup requires alignment step. Please enable alignment step using --runSTAR!"
      }

      // Include required modules.
      include MarkDup from '../NextflowModules/Sambamba/0.7.1/MarkDup.nf' params(mem:params.sambambamarkdup.mem, optional:params.markdup.toolOptions)

      // Run workflow.
      workflow markdup_mapping {
          take:
            bams_in
          main:
            /* Run mapping on a per sample per lane basis */
            MarkDup(bams_in)
          emit:
            bams = MarkDup.out
      }

      markdup_mapping(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, sample_id, bams, bai] })

    }

    // Generate post-QC plots.
    if ( params.runPostQC) {
      // Input validation.
      if (!params.runSTAR){
        exit 1, "runPostQC requires alignment step. Please enable alignment step using --runSTAR!"
      }

      // Check / generate BED files needed for TIN.
      if (params.genome_bed) {
        // Load supplied BED12 file.
        genome_bed = Channel
              .fromPath(params.genome_bed, checkIfExists: true)
              .ifEmpty { exit 1, "Supplied BED12 file not found: ${params.genome_bed} (--genome_bed)"}
      } else if ( !params.genome_bed) {

        // Include required modules.
        include GtfToGenePred from './NextflowModules/ucsc/377/gtfToGenePred/GtfToGenePred.nf' params(params)
        include GenePredToBed from './NextflowModules/ucsc/377/genePredToBed/GenePredToBed.nf' params(params)

        // Generate BED12 file.
        GtfToGenePred ( genome_gtf)
        GenePredToBed ( GtfToGenePred.out )
        genome_bed = GenePredToBed.out
      }

      // Include required modules.
      include post_mapping_QC from './sub-workflows/post_mapping_QC.nf' params(params)

      // Run workflow.
      post_mapping_QC(mapped.map { sample_id, bams, unmapped, log1, log2, tab, bai -> [sample_id, bams, bai] }, genome_bed.collect())

      // Logging
      logs_PostQC =  post_mapping_QC.out[1].map { it[1] }.mix(post_mapping_QC.out[0].map { it[1] })

    }else{

      // Create empty log.
      logs_PostQC = Channel.empty()

    }


    // Read counting - FeatureCounts / SubRead
    if ( params.runFeatureCounts) {

      // Include required modules.
      include FeatureCounts from './NextflowModules/subread/2.0.0/FeatureCounts.nf' params(
        optional:params.fc.toolOptions,
        singleEnd: params.singleEnd,
        stranded: params.stranded,
        unstranded: params.unstranded,
        revstranded: params.revstranded,
        fc_group_features: params.fc_group_features,
        fc_count_type: params.fc_count_type,
        fc_group_features_type: params.fc_group_features_type,
        fc_extra_attributes : params.fc_extra_attributes,
        gencode: params.gencode
      )

      // Run workflow.
      FeatureCounts(run_name, AlignReads.out.map { it[1] }.collect(), genome_gtf.collect())

      if ( params.normalize_counts ) {
        // Include required modules.
        include EdgerNormalize as fc_norm from './utils/bioconductor/edger/3.28.0/normalize.nf' params(
          tool:"fc"
        )

        fc_norm( run_name, FeatureCounts.out.map { it[1] } )

      }

      // Logging
      logs_FeatureCounts = FeatureCounts.out.map { it[3]}

    }else{

      // Create empty log.
      logs_FeatureCounts = Channel.empty()

    }

    // Aligment and quantification using Salmon.
    if ( params.runSalmon ) {

      // Input validation.
      if (!params.transcripts_fasta) {
        exit 1, "Transcript fasta not found for Salmon. Please provide the correct path! (--transcripts_fasta)"
      }else{
        transcripts_fasta = Channel
            .fromPath(params.transcripts_fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Supplied transcript FASTA file not found: ${params.transcripts_fasta} (--transcripts_fasta)"}
      }

      // Load / generate Salmon index.
      if ( params.salmon_index ) {
         salmon_index = Channel
              .fromPath(params.salmon_index, checkIfExists: true)
              .ifEmpty { exit 1, "Supplied Salmon index not found: ${params.salmon_index} (--salmon_index)"}
      } else {

        // Include required modules.
        include Index as SalmonIndex from './NextflowModules/Salmon/1.2.1/Index.nf' params(
          gencode: params.gencode
        )

        // Generate Salmon Index.
        SalmonIndex ( transcripts_fasta )
        salmon_index = SalmonIndex.out
      }

      // Include required modules.
      include QuantMerge from './NextflowModules/Salmon/1.2.1/QuantMerge.nf' params(
        optional: params.salmon_quantmerge.toolOptions
      )

      include Quant from './NextflowModules/Salmon/1.2.1/Quant.nf' params(
        singleEnd: params.singleEnd,
        stranded: params.stranded,
        unstranded: params.unstranded,
        revstranded: params.revstranded,
        saveUnaligned: params.saveUnaligned,
        optional: params.salmon_quant.toolOptions
      )

      include mergeFastqLanes from './NextflowModules/Utils/mergeFastqLanes.nf' params(params)

      // Run workflow.
      Quant ( mergeFastqLanes (fastqs_transformed ), salmon_index.collect() )
      QuantMerge ( Quant.out.map { it[1] }.collect(), run_name )

      // Logging
      logs_Salmon = Quant.out.map { it[1] }

    }else{

      // Create empty log.
      logs_Salmon = Channel.empty()

    }

    // ToDo - Split from GATK workflow.
    // Perform BQRS - GATK
    if( params.runBQRS) {

      // Include required modules.

    }

    // GATK (v4) - Germline calling
    if ( params.runGermlineCallingGATK ) {

      // Input validation
      if (!params.genome_dict) {
          exit 1, "Genome dictionary not found for GATK, please provide the correct path! (--genome_dict)"
      }

      if (!params.genome_index) {
          exit 1, "Genome index not found for GATK, please provide the correct path! (--genome_index)"
      }

      if(!params.runSTAR && params.runMarkDup){
        exit 1, "--runGermlineCallingGATK required aligned and duplicates marked BAM files! (--runSTAR & --runMarkDup)"
      }

      // Load / Generate Scatter interval list.
      if (params.scatter_interval_list ) {
        scatter_interval_list = Channel
          .fromPath( params.scatter_interval_list, checkIfExists: true)
          .ifEmpty { exit 1, "Supplied scatter intervals list not found: ${params.scatter_interval_list} (--scatter_interval_list)"}
      } else if ( !params.scatter_interval_list ) {

        // Include required modules.
        include CreateIntervalList from './NextflowModules/Utils/CreateIntervaList.nf' params(params)

        // Generate scatter interval list.
        genome_dict = Channel
              .fromPath( params.genome_dict, checkIfExists: true)
              .ifEmpty { exit 1, "Supplied genome dictionary not found: ${params.genome_dict} (--genome_dict)"}
        genome_index = Channel
              .fromPath(params.genome_fasta + '.fai', checkIfExists: true)
              .ifEmpty { exit 1, "Supplied FASTA does not has accompanying FAI file in the same folder: ${params.genome_fasta}.fai (--genome_fasta)"}

        // Run workflow.
        CreateIntervalList( genome_index, genome_dict )
        scatter_interval_list = CreateIntervalList.out

      }

      // Include required modules.
      include SplitIntervals from './NextflowModules/GATK/4.1.3.0/SplitIntervals.nf' params(
        optional: params.splitintervals.toolOptions
      )
      include gatk4_bqsr from './sub-workflows/gatk4_bqsr.nf' params(params)
      include gatk4_hc from './sub-workflows/gatk4_hc.nf' params(params)
      include SplitNCigarReads from './NextflowModules/GATK/4.1.3.0/SplitNCigarReads.nf' params(genome_fasta:params.genome_fasta)

      // Run workflow.
      SplitIntervals( 'no-break', scatter_interval_list)
      SplitNCigarReads(markdup_mapping.out)

      if ( params.runGATK4_BQSR) {

        //Perform BSQR
        gatk4_bqsr(SplitNCigarReads.out, SplitIntervals.out.flatten())
        gatk4_hc(gatk4_bqsr.out[0], SplitIntervals.out.flatten(), run_name)

      } else {
          gatk4_hc(SplitNCigarReads.out, SplitIntervals.out.flatten(), run_name)
      }
    }

    // Combine reports - MultiQC.
    if ( params.runMultiQC ) {

      // Include required modules.
      include multiqc_report from './sub-workflows/multiqc_report.nf' params(params)

      // Run workflow.
      multiqc_report( logs_FastQC,
                      logs_TrimGalore,
                      logs_SortMeRNA,
                      logs_STAR,
                      logs_PostQC,
                      logs_FeatureCounts,
                      logs_Salmon )
    }

}

// End of the workflow.
