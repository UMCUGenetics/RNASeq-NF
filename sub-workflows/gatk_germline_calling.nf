if (!params.genome_fasta) {
  exit 1, "GATK requires a genome fasta file. Please provide the correct filepath! (--genome_fasta)"
} 
//Import GATK modules
include { HaplotypeCaller } from params.nextflowmodules_path+'/GATK/4.3.0.0/HaplotypeCaller.nf' params ( mem:params.haplotypecaller.mem,
                                                                                           compress: params.compress,
                                                                                           genome_fasta:params.genome_fasta,
                                                                                           optional: params.options.GATK4_HaplotypeCaller )
include { GenotypeGVCFs } from params.nextflowmodules_path+'/GATK/4.3.0.0/GenotypeGVCFs.nf' params(mem: "", 
                                                                                            genome_fasta: "${params.genome_fasta}", 
                                                                                            genome_dbsnp: "${params.genome_dbsnp}", 
//                                                                                            genome_dbsnp: "", 
                                                                                            compress: "${params.compress}")
include { CombineGVCFs } from params.nextflowmodules_path+'/GATK/4.3.0.0/CombineGVCFs.nf' params(mem: "", 
                                                                                            genome_fasta: "${params.genome_fasta}", 
                                                                                            compress: "${params.compress}")
include { VariantFiltration } from params.nextflowmodules_path+'/GATK/4.3.0.0/VariantFiltration.nf' params( mem:params.variantfiltration.mem,
                                                                                              genome_fasta:params.genome_fasta,
                                                                                              gatk_snp_filter: params.gatk_snp_filter, 
                                                                                              gatk_indel_filter: params.gatk_indel_filter, 
                                                                                              gatk_rna_filter: params.gatk_rna_filter,
                                                                                              compress: params.compress,
                                                                                              optional: params.options.GATK4_VariantFiltration )
include { MergeVCFs as MergeVCF } from params.nextflowmodules_path+'/GATK/4.3.0.0/MergeVCFs.nf' params( mem:params.mergevcf.mem,
                                                                                                        compress: params.compress)
include { BaseRecalibrationTable } from params.nextflowmodules_path+'/GATK/4.3.0.0/BaseRecalibrationTable.nf' params( mem:params.baserecalibrator.mem,
                                                                                                        optional:params.options.GATK4_BQRS,
                                                                                                        genome_known_sites:params.genome_known_sites,
                                                                                                        genome_fasta:params.genome_fasta )
include { GatherBaseRecalibrationTables } from params.nextflowmodules_path+'/GATK/4.3.0.0/GatherBaseRecalibrationTables.nf' params( mem:params.gatherbaserecalibrator.mem )
include { BaseRecalibration } from params.nextflowmodules_path+'/GATK/4.3.0.0/BaseRecalibration.nf' params ( mem:params.applybqsr.mem,
                                                                                               genome_fasta:params.genome_fasta )
include { MergeBams } from params.nextflowmodules_path+'/Sambamba/0.8.2/MergeBams.nf' params(mem: "${params.mergebams.mem}")
include { SplitIntervals } from params.nextflowmodules_path+'/GATK/4.3.0.0/SplitIntervals.nf' params( optional: params.options.GATK4_SplitIntervals )
include { SplitNCigarReads } from params.nextflowmodules_path+'/GATK/4.3.0.0/SplitNCigarReads.nf' params( genome_fasta:params.genome_fasta)                     
include { CreateIntervalList } from params.nextflowmodules_path+'/Utils/CreateIntervaList.nf' params( params )
//Sambamba modules
include { Markdup } from params.nextflowmodules_path+'/Sambamba/0.8.2/Markdup.nf' params( mem:params.sambambamarkdup.mem )
include { Flagstat as Flagstat_markdup } from params.nextflowmodules_path+'/Sambamba/0.8.2/Flagstat.nf' params( params )


workflow gatk_germline_calling {
    take:
      run_id
      bam_file

    main:
        //Check for Scatter intervallist
        if ( params.scatter_interval_list ) {
            scatter_interval_list = Channel
              .fromPath( params.scatter_interval_list, checkIfExists: true)
              .ifEmpty { exit 1, "Genome scatter intervals not found: ${params.scatter_interval_list}"}
        } else if ( !params.scatter_interval_list ) {
            genome_index = Channel
                .fromPath(params.genome_fasta + '.fai', checkIfExists: true)
                .ifEmpty { exit 1, "Genome fai file not found: ${params.genome_fasta}.fai"}
            genome_dict = Channel
                .fromPath( params.genome_dict, checkIfExists: true)
                .ifEmpty { exit 1, "Genome dict not found: ${params.genome_dict}"}
            CreateIntervalList( genome_index, genome_dict )
            scatter_interval_list = CreateIntervalList.out.genome_interval_list
        }
        //Scatter intervals
        SplitIntervals( 'no-break', scatter_interval_list)
        scatter_intervals = SplitIntervals.out.flatten()
        //Sambamba Markdup
        Markdup( bam_file )
        Flagstat_markdup( Markdup.out.map {sample_id, rg_id, bam, bai -> 
			                  [sample_id, bam, bai]} )
        //NCigar split
        SplitNCigarReads( Markdup.out.map {sample_id, rg_id, bam, bai ->
                                          [sample_id, bam, bai]} )
        final_bam = SplitNCigarReads.out.bam_file
        //Perform BQSR
        if ( params.runGATK4_BQSR  ) {
            if ( !params.genome_known_sites ) {
                 exit 1, "BQSR requires known variant sites for recalibration (--genome_known_sites). "
            }
            BaseRecalibrationTable(final_bam.combine(scatter_intervals))
            GatherBaseRecalibrationTables(BaseRecalibrationTable.out.groupTuple())
            //Perform BQSR
            BaseRecalibration(
              final_bam
                .combine(GatherBaseRecalibrationTables.out, by:0)
                .combine(scatter_intervals)
            )
            //Merge recalibrated bams
              MergeBams(
              BaseRecalibration.out
                .groupTuple()
                .map{ [it[0],it[2],it[3]] }
            )
            //Set final bam file to recalibrated bam
            final_bam = MergeBams.out
        }
        //      
        HaplotypeCaller(final_bam.combine( scatter_intervals) )

        //Legacy mode generates a single VCF per sample
        if( params.legacyMode ){
            //Merge scattered vcf chunks/sample
            MergeVCF(
              HaplotypeCaller.out.groupTuple(by:[0]).map{
                sample_id, intervals, vcfs, idxs, interval_files ->
                [sample_id, vcfs, idxs] }
            )
        } else {
            //By default now one VCF is generated for the run
            gvcf_per_interval = HaplotypeCaller.out
                .groupTuple(by:[1])
                .map{ sample_ids, interval, gvcfs, idxs, interval_files ->
                      [run_id, interval, gvcfs, idxs, interval_files[0]]
                }
        
            //Combine GVCFs per interval (all samples per interval)
            CombineGVCFs( gvcf_per_interval )
            GenotypeGVCFs(CombineGVCFs.out)

            //merge raw/unfiltered vcfs
            MergeVCF(
                GenotypeGVCFs.out.groupTuple().map{ 
                    run_id, interval, gvcfs, idxs, intervalfile ->
                    [run_id, gvcfs, idxs ] }
            )
        }

        //Filter raw vcf files/sample
        VariantFiltration( MergeVCF.out.map{
          sample_id, vcfs, idxs -> [sample_id, run_id, "RNA", vcfs, idxs] }
    )
      
    emit:
      bam_recal = MergeBams.out
      bam_markdup = Markdup.out
      markdup_flagstat = Flagstat_markdup.out
      bam_ncigar =  SplitNCigarReads.out	
      vcf_filter = VariantFiltration.out
      bqsr_table = GatherBaseRecalibrationTables.out

}
