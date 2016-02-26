---
title: "Alignment-free expression estimation using Sailfish"
author: "Radhika Khetani"
date: ""
---

Contributors: Radhika Khetani, Meeta Mistry

Approximate time: 3 hours

## Learning Objectives

* Understand the concept of "Pseudocounts"
* Understand how to use Sailfish to generate Psuedocounts
* Learn how to perform differential gene expression on Psuedocounts

## What is Sailfish



## Analysis steps for Sailfish

## Getting the abundance estimates

    mkdir ~/ngs_course/rnaseq/sailfish
    cd ~/ngs_course/rnaseq/sailfish
    
    # sailfish index sailfish index -p <num of cores> -k <kmer size> -t <fasta of gene sequences> -o <folder name>
    
    export PATH=/groups/bcbio/bcbio/anaconda/bin:/opt/bcbio/local/bin:$PATH
    
    sailfish quant -i /groups/hbctraining/sailfish-run/sailfish.ensembl2.idx/ \
    -l SR \
    -r ngs_course/rnaseq/data/untrimmed_fastq/Mov10_oe_1.subset.fq \
    --useVBOpt \
    -o Mov10_oe_1.subset.sailfish
    
    less Mov10_oe_1.subset.sailfish/quant.sf
  
## Converting to psuedo-counts

## Using DESeq2 for DGE analysis  
  
  
***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
