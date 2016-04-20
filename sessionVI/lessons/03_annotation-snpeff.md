---
title: "Variant annotation with snpEff"
author: "Meeta Mistry, Mary Piper"
date: "Thursday, April 28th, 2016"
---

Approximate time: 90 minutes

## Learning Objectives:

* ()
* ()



## Annotating variants

Variant annotation is a crucial step in linking sequence variants with changes in phenotype. Annotation results can have a strong influence on the ultimate conclusions of disease studies. Incorrect or incomplete annotations can cause researchers both to overlook potentially disease-relevant DNA variants and to dilute interesting variants in a pool of false positives. At this stage, we have a large tab-delimited file containing loci at which a variation was found in the sample DNA sequence relative to the reference. We have filtered out these variations (also referred to as 'variant calls') to keep only those we are highly confident in, and now need to find out more. We can do this by comparing our variants against known variants, and also use genome annoatations to help predict information about our variants. 

![var_calling_workflow](../img/variant_calling_workflow.png)

### Setting up

For this section we are going to need to copy over some reference data required for annotation. First, move into `var-calling` directory and then copy over the required data. *NOTE: one is a directory and the other is a zipped file*

```
$ cd ngs_course/var-calling
$ cp -r /groups/hbctraining/ngs-data-analysis2016/var-calling/reference_data/snpeff data/reference_data
$ cp /groups/hbctraining/ngs-data-analysis2016/var-calling/reference_data/dbsnp.138.chr20.vcf.gz data/reference_data/

```

Let's also create a new directory for the annotation steps:

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



## Comparing against known variants 

Variant annotation is the process of assigning information to DNA variants. There are many different types of information that could be associated with variants, and a first commonly use resource is to compare against variants that have previously been described. [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) 

## Functional annotation with SnpEff

One fundamental level of variant annotation involves categorising each variant based on its relationship to coding sequences in the genome and how it may change the coding sequence and affect the gene product. To do this we will be using a tool called [SnpEff](http://snpeff.sourceforge.net/), a variant effect predictor program. 

Our understanding of the protein-coding sequences in the genome is summarised in the set of transcripts we believe to exist. Thus, variant annotation depends on the set of transcripts used as the basis for annotation. The widely used annotation databases and browsers – ENSEMBL, RefSeq, UCSC – contain sets of transcripts that can be used for variant annotation, as well as a wealth of information of many other kinds as well, such as ENCODE data about the function of non-coding regions of the genome. 

Thus, the first step to using SnpEff will involve building a database based on annotations retrieved from a public database. **We have generated this database for you, using Ensembl annotations provided within the snpEff framework.** You will need to copy over this database directory into your `reference` directory:


Now we can run the snpEff command 





***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*




