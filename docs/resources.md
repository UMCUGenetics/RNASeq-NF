### 1. Reference sequence and annotation (required)
Download your prefered reference genome from [Ensembl](https://www.ensembl.org/index.html), [Illumina iGenomes](https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html) etc. and stored them in a dedicated resource folder. This should include at least the reference sequence (.fasta) and genome annotation (.gtf). 

### 2. Reference transcript sequence (optional)
In case you intend to run transcript quantification with Salmon, download the cDNA sequences of interest from [Ensembl](https://www.ensembl.org/index.html). Alternatively, the transcript fasta can be generated with [RSEM](http://deweylab.biostat.wisc.edu/rsem/rsem-prepare-reference.html) directly from genomic .fasta or .gtf inputs. 

```
rsem-prepare-reference --gtf Homo_sapiens.GRCh37.74.gtf \
                       --star \ 
                        Homo_sapiens.GRCh37.fasta /Homo_sapiens_GRCh37.transcripts
```

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

These files are not required if BQSR recalibration is disabled.
