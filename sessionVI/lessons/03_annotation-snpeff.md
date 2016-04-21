---
title: "Variant annotation with snpEff"
author: "Meeta Mistry, Mary Piper"
date: "Thursday, April 28th, 2016"
---

Approximate time: 90 minutes

## Learning Objectives:

* Adding infomration on known SNPs to our VCF
* Adding functional information to variants in the VCF
* Understanding where the annotation is added in the VCF format



## Annotating variants

Variant annotation is a crucial step in linking sequence variants with changes in phenotype. Annotation results can have a strong influence on the ultimate conclusions of disease studies. Incorrect or incomplete annotations can cause researchers both to overlook potentially disease-relevant DNA variants and to dilute interesting variants in a pool of false positives. 

![var_calling_workflow](../img/variant_calling_workflow.png)

At this stage, we have a large tab-delimited file containing loci at which a variation was found in the sample DNA sequence relative to the reference. We have filtered out these variations (also referred to as 'variant calls') to keep only those we are highly confident in, and now need to find out more. We can do this by **comparing our variants against known variants, and also use genome annoatations to help predict information about our variants.** 

<img src="../img/prioritize.png" width="300">



### Setting up

For this section we are going to need to copy over some reference data required for annotation. First, move into `var-calling` directory and then copy over the required data.

```
$ cd ngs_course/var-calling
$ cp /groups/hbctraining/ngs-data-analysis2016/var-calling/reference_data/dbsnp.138.chr20.vcf.gz* data/reference_data/

```

Let's also create a new directory for the results of our annotation steps:

```
$ mkdir results/annotation

```
Your directory structure should now look something like this:

```
~/
├── ngs_course/
    ├── var-calling/
        ├── data/
            ├── untrimmed_fastq/
            ├── reference_data
        	    ├── chr20.fa
        	    ├── dbsnp.138.chr20.vcf.gz
                └── snpeff
        ├── results/
     		├── bwa
     		├── variants
            └── annotation

```



## Annotation with known variants 

Variant annotation is the process of assigning information to DNA variants. There are many different types of information that can be associated with variants, and a first commonly used resource is using databases which contain variants that have previously been described. One popular example is [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/),a free, public archive for genetic variation within and across different species. It is hosted by NCBI in collaboration with NHGRI and although the name implies SNPs; it actually includes range of molecular variation.

<img src="../img/dbsnp.png" width="300">

To add dbSNP information you need to download the organism specific data using their [FTP download](ftp://ftp.ncbi.nih.gov/snp/) site. **We have already done this for you** and was the zip file that you copied over into your `reference_data` folder. 

To annotate our data with dbSNP information we wil be using [`bcftools`](https://samtools.github.io/bcftools/), a command-line utility for variant calling and manipulating VCF files and its binary counterpart BCF files. It is a part of the `samtools` project, a tool that we are by now pretty familiar with. 

The `bcftools annotate` command allows the user to add or remove annotations. The annotation we wish to add must be a Bgzip-compressed and tabix-indexed file (usually VCF or BED format), as should the file that we are annotating:

```
$ bgzip results/variants/na12878_q20.recode.vcf 
$ tabix results/variants/na12878_q20.recode.vcf.gz
```

When running `bcftools annotate`, we also need to specify the column(s) to carry over from the annotation file, which in our case is ID.

```
$ bcftools annotate -c ID -a data/reference_data/dbsnp.138.chr20.vcf.gz results/variants/na12878_q20.recode.vcf.gz > results/annotation/na12878_q20_annot.vcf
```
Take a quick peek at the new VCF file that was generated using `less`. You should now see in the ID column `rs` ids which correspond to identifiers from dbSNP. For the variants  that are not known, you will find the `.` in place of an ID indicating novelty.

 

## Functional annotation with SnpEff

One fundamental level of variant annotation involves categorising each variant based on its relationship to coding sequences in the genome and how it may change the coding sequence and affect the gene product. To do this we will be using a tool called [SnpEff](http://snpeff.sourceforge.net/), a variant effect predictor program. 

Our understanding of the protein-coding sequences in the genome is summarised in the set of transcripts we believe to exist. Thus, variant annotation depends on the set of transcripts used as the basis for annotation. The widely used annotation databases and browsers – ENSEMBL, RefSeq, UCSC – contain sets of transcripts that can be used for variant annotation, as well as a wealth of information of many other kinds as well, such as ENCODE data about the function of non-coding regions of the genome. 

Thus, the first step to using SnpEff will involve building a database based on annotations retrieved from a public database. **We have generated this database for you, using Ensembl annotations provided within the snpEff framework.** You will need to copy over this database directory into your `reference` directory:


Now we can run the snpEff command 





***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*




