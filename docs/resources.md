### 1. Reference sequence and annotation (required)
Download your prefered reference genome from [Ensembl](https://www.ensembl.org/index.html), [Illumina iGenomes](https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html) etc. and stored them in a dedicated resource folder. This should include at least the reference sequence (.fasta) and genome annotation (.gtf). 

### 2. Reference transcript sequence (optional)
For transcript quantificaton with Salmon, download the cDNA sequences of interest from [Ensembl](https://www.ensembl.org/index.html).

For example;          
```
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz
```
The fasta will be used to build the Salmon transcriptome index for alignment-free quantification.


### 3. GATK resource bundle (optional)  

Download the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle) for your genome build of interest if you intend to perform SNP/Indel calling on RNASeq data. For example, the GRCh37 resource files can be obtained as shown below.

```
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.idx.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.idx.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.idx.gz
```

### 4. SortMeRNA ribosomal reference files (optional)
The locations of the rRNA databases are stored in the /resources/sortmerna-db-default.txt file. The fasta files are downloaded from GitHub and staged in the pipeline working directory. Use the offline version in case of connection problems or unavailability. 
