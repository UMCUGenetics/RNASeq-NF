include { MarkdownToPdf } from '../NextflowModules/Pandocker/21.02/MarkdownToPdf.nf' params( params )
include { FillMdTemplate } from '../utils/FillRNASeqMdTemplate.nf' params( params )    
include { CreateTmpValueMap } from '../utils/FillRNASeqMdTemplate.nf' params( params )    

workflow generate_rnaseq_report {
    take:
        run_name
        md_template
        value_map_file

    main:
        if ( !value_map_file ){
            CreateTmpValueMap( run_name )
            value_map_file = CreateTmpValueMap.out.value_map_file
        }
        FillMdTemplate( run_name, md_template, value_map_file )
        report_md_file = FillMdTemplate.out.md_file

        MarkdownToPdf( report_md_file )
      
    emit:
        pdf_report =  MarkdownToPdf.out
        md_report = report_md_file
}
