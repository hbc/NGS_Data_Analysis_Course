---
title: "Counting reads"
author: "Meeta Mistry, Bob Freeman, Radhika Khetani"
date: "Friday, February 12, 2016"
---

Approximate time: 

## Learning Objectives:

* learn about tools that generate read counts, a measure of gene expression
* learn how to use the featureCounts tool to generate a count matrix for statistical analyses

## Counting reads as a measure of gene expression
<img src="../img/rnaseq_workflow.png" width="400">

Once we have our reads aligned to the genome, the next step is to count how many reads have mapped to each gene. There are many tools that can use BAM files as input and output the number of reads (counts) associated with each feature of interest (genes, exons, transcripts, etc.). There are 2 commonly used counting tools, [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) and [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html). 

* The above tools only report the "raw" counts of reads that map to a single location (uniquely mapping) and are best at counting at the gene level. Essentially, total read count associated with a gene (*attribute*) = the sum of reads associated with each of the exons (*feature*) that "belong" to that gene.

* There are other tools available that are able to account for multiple transcripts for a given gene. In this case the counts are not whole numbers, but have fractions. In the simplest example case, if 1 read is associated with 2 transcripts, it can get counted as 0.5 and 0.5 and the resulting count for that transcript is not a whole number.

> **Input for counting**: BAM files + GTF file.
> Simply speaking, the genomic coordinates of where the read is mapped (BAM) are cross-referenced with the genomic coordinates of whichever feature you are interested in counting expression of (GTF), it can be exons, genes or transcripts.

<img src="../img/count-fig1.png" width="600">

> **Output of counting**: A count matrix, with genes as rows and samples are columns. 
> These are the "raw" counts and will be used in statistical programs downstream for differential gene expression.

<img src="../img/count-matrix.png" width="500">

### Counting using featureCounts
Today, we will be using the featureCounts tool to get the *gene* counts, since this tool is accurate and it is relatively easy to use. This tool only counts reads that are mapping to a single location (uniquely mapping) and follows the scheme in the figure below for assigning reads to a gene/exon. 

<img src="../img/union.png" width="400">

(figure adapted from http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)

featureCounts can take into account whether your data are stranded or not. If strandedness is specified then, in addition to considering the genomic coordinates, it will also take the strand into account for counting.

First things first, start an interactive session with 4 cores
	
	$ bsub -Is -n 4 -q interactive bash

Now, change directories to your rnaseq directory and start by creating 2 directories, (1) a directory for the output and (2) a directory for just the bam files we generated yesterday:

	$ cd ~/ngs_course/rnaseq/
	$ mkdir results/counts results/STAR/bam

featureCounts is not available as a module on Orchestra, but we can add the path for it to our PATH variable. 

	$ export PATH=/opt/bcbio/local/bin:$PATH

> Remember that this export command is only valid for this interactive session. If you want to make sure that the tool is available to you all the time, add the above command to your `~/.bashrc` or your `~/.bash_profile` files.

What options/parameters are available to us for this tool?

	```$ featureCounts 

	  Version 1.4.4

	  Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 
	  Required arguments:

	  -a <string>         Name of an annotation file. GTF format by default. See -F 
	                      option for more formats.

	  -o <string>         Name of the output file including read counts. A separate 
	                      file including summary statistics of counting results is 
	                      also included in the output (`<string>.summary')

	  input_files         List of input files in BAM or SAM format. Users do not 
	                      need to specify it is BAM or SAM.
	  
	  Optional arguments:

	  .
	  .
	  .
	  .```

Now

	featureCounts -T 4 -a ~/ngs_course/unix_lesson/reference_data/chr1-hg19_genes.gtf \
		-o ~/ngs_course/unix_lesson/rnaseq/results/counts/counts.txt \
		-s 2 \
		~/ngs_course/unix_lesson/rnaseq/results/STAR/*bam


#### Exercise
Take a look at the end of the file using the `tail` command. You should see a summary of how the reads were classified. 

1. How many reads were assigned as no_feature? Why would they be classified this way?
2. How many reads were found to map to multiple locations?






