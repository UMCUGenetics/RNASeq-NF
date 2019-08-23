#!/usr/bin/env nextflow

params.run_fastqc = true
params.run_star = true

nextflow.preview.dsl=2


if (params.endness == "paired"){
  //Return paired-end channel
  Channel.fromPath( file(params.samplesheet) )
    .splitCsv(header: true, sep: '\t')
    .map{row ->
        def sample = row['sampleid']
        def reads1 = row['R1'].tokenize( ',' ).collect { file(it) }
        def reads2 = row['R2'].tokenize( ',' ).collect { file(it) }
        return [ sample, reads1, reads2 ]
    }
    .tap{samples_R1_R2_fastq} // set of all fastq R1 R2 per sample
    .map { sample, reads1, reads2 ->
        return [ sample, [reads1, reads2 ].flatten() ]
    }
    .tap{samples_all_fastq} // set of all fastq per sample
}
//Return single ended channel
else {
   Channel.fromPath( file(params.samplesheet) )
    .splitCsv(header: true, sep: '\t')
    .map{row ->
        def sample = row['sampleid']
        def reads1 = row['R1'].tokenize( ',' ).collect { file(it) }
        return [ sample, reads1  ]
    }
    .tap{samples_all_fastq; samples_all_fastq_star } // set of all fastq per sam
}

if (params.run_fastqc){
    include FastQC from 'NextflowModules/FastQC/0.11.8/FastQC.nf' params(params)
    FastQC(samples_all_fastq)
}

if (params.run_star) {
    index = Channel.fromPath(params.star_index)
    include AlignReads from 'NextflowModules/STAR/2.4.2a/AlignReads.nf' params(params)
    include featureCount from 'NextflowModules/HTseq/featureCount.nf' params(params)
    aligned = AlignReads(samples_R1_R2_fastq, index.collect())
    featureCount(aligned[0])
}












