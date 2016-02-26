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

[Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/index.html) and it's more recent "upgrade" [Salmon](https://combine-lab.github.io/salmon/), are based on the philosophy of lightweight algorithms. They use the sequence of genes or transcripts as input, and do not align the whole read. Instead it's a 2-step process:

a. they first evaluate the sequences for all possible unique sequences of lenght k (kmer) in the transcriptome (genes/transcripts).

b. then they count the number of times those kmers appear in the sequenced data, i.e. the fastq. This count information is then used to eestimate the abundance of each gene or transcript. 

<img src="../img/nbt.2862-F1.jpg" width="400">

## Analysis steps for Sailfish

    mkdir ~/ngs_course/rnaseq/sailfish
    cd ~/ngs_course/rnaseq/sailfish
    export PATH=/groups/bcbio/bcbio/anaconda/bin:/opt/bcbio/local/bin:$PATH
    
As you can imagine from the above schematic, there are 2 steps in the analysis too:
1. "Index" the transcriptome (transcripts or genes) as follows:
    
    # sailfish index sailfish index -p <num of cores> -k <kmer size> -t <fasta of gene sequences> -o <folder name>
> We are not going to run this in class, but it only takes a few minutes.

2. Get the abundance using the quantification step as follows:

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
