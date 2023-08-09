include { MarkdownToPdf } from params.nextflowmodules_path+'/Pandocker/21.02/MarkdownToPdf.nf' params( params )
include { MarkdownToHtml } from params.nextflowmodules_path+'/Pandocker/21.02/MarkdownToHtml.nf' params( params )
include { FillMdTemplate } from '../utils/FillRNASeqMdTemplate.nf' params( params )    

workflow generate_rnaseq_report {
    take:
        run_name
        md_template

    main:
        report_name=run_name+"_RNASeq_report"
        FillMdTemplate( report_name, md_template )
        report_md_file = FillMdTemplate.out.md_file

        MarkdownToPdf( report_md_file )
        MarkdownToHtml( report_md_file )
      
    emit:
        pdf_report =  MarkdownToPdf.out
        md_report = report_md_file
}
