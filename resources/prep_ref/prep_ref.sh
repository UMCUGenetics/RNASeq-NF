#STAR index
STAR --runMode genomeGenerate --genomeDir <out_dir> \
        --genomeFastaFiles <fasta> \
        --runThreadN 4 \
        --sjdbGTFfile <gtf>

#Calculate exon gene lengths
Rscript --vanilla exon_lengths.R <gtf> <genome>_exon_gene_sizes.txt

#GTF to Bed12 conversion (https://hgdownload.cse.ucsc.edu/admin/exe)
gtfToGenePred <gtf> out.genePhred
genePredToBed out.genePhred out.bed12
sort -k1,1 -k2,2n out.bed12 > out.bed12.sorted.bed

#Salmon transcripts for quantification
$RSEM_ROOT/rsem-prepare-reference --gtf <gtf> --p 8 --star <fasta> <rsem_out>


