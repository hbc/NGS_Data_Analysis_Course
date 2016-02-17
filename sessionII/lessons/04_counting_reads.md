---
title: "Counting reads"
author: "Meeta Mistry, Bob Freeman, Radhika Khetani"
date: "Friday, February 12, 2016"
---

Approximate time: 

## Learning Objectives:

* learn how to use the featureCounts tool to generate a count matrix for statistical analyses

## Counting reads as a measure of gene expression
<img src="../img/rnaseq_workflow.png" width="400">

Once we have our reads aligned to the genome, the next step is to count how many reads have mapped to each gene. There are many tools that can use BAM files as input and output the number of reads (counts) associated with each feature of interest (genes, exons, transcripts, etc.). There are 2 commonly used counting tools, [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) and [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html). 

* The above tools only report the "raw" counts of reads that map to a single location (uniquely mapping) and are best at counting at the gene level. Essentially, total read count associated with a gene (*meta-feature*) = the sum of reads associated with each of the exons (*feature*) that "belong" to that gene.

* There are other tools available that are able to account for multiple transcripts for a given gene. In this case the counts are not whole numbers, but have fractions. In the simplest example case, if 1 read is associated with 2 transcripts, it can get counted as 0.5 and 0.5 and the resulting count for that transcript is not a whole number.

**Input for counting**: BAM files + GTF file.
Simply speaking, the genomic coordinates of where the read is mapped (BAM) are cross-referenced with the genomic coordinates of whichever feature you are interested in counting expression of (GTF), it can be exons, genes or transcripts.

<img src="../img/count-fig1.png" width="600">

**Output of counting**: A count matrix, with genes as rows and samples are columns. 
These are the "raw" counts and will be used in statistical programs downstream for differential gene expression.

<img src="../img/count-matrix.png" width="300">

### Counting using featureCounts
Today, we will be using the featureCounts tool to get the *gene* counts, since this tool is accurate and it is relatively easy to use. This tool only counts reads that are mapping to a single location (uniquely mapping) and follows the scheme in the figure below for assigning reads to a gene/exon. 

<img src="../img/union.png" width="300">

featureCounts can also take into account whether your data are **stranded** or not. If strandedness is specified then, in addition to considering the genomic coordinates, it will also take the strand into account for counting. If your data are stranded, please do specify it.

First things first, start an interactive session with 4 cores
	
	$ bsub -Is -n 4 -q interactive bash

Now, change directories to your rnaseq directory and start by creating 2 directories, (1) a directory for the output and (2) a directory for just the bam files we generated yesterday:

	$ cd ~/ngs_course/rnaseq/
	$ mkdir results/counts results/STAR/bams
	
Let's move over the bam files over to the `results/STAR/bams` directory
	
	$ mv ~/ngs_course/rnaseq/results/STAR/*fq_Aligned*bam ~/ngs_course/rnaseq/results/STAR/bams

featureCounts is not available as a module on Orchestra, but we can add the path for it (parent directory) to our PATH variable. 

	$ export PATH=/opt/bcbio/local/bin:$PATH

> Remember that this export command is only valid for this interactive session. If you want to make sure that the tool is available to you all the time, add the above command to your `~/.bashrc` or your `~/.bash_profile` files.

How do we use this tool, what is the command and what options/parameters are available to us?

	$ featureCounts

So, it looks like the usage is `featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... `, where `-a`, `-o` and input files are required.

We are going to use the following options:

`-T 4` # specify 4 cores
`-s 2` # these data are "reverse"ly stranded

and the following are the values for the required parameters:

`-a ~/ngs_course/rnaseq/data/reference_data/chr1-hg19_genes.gtf` # required option. Specify path to GTF
`-o ~/ngs_course/rnaseq/results/counts/Mov10_featurecounts.txt` #  required option. Specify path to, and name of the text output (count matrix)
`~/ngs_course/rnaseq/results/STAR/bams/*bam` # the list of all the bam files we want to collect count information for

We are now going to run this in the interactive session using 4 cores.

	$ featurecounts -T 4 -s 2\ 
	  -a ~/ngs_course/unix_lesson/reference_data/chr1-hg19_genes.gtf \
	  -o ~/ngs_course/unix_lesson/rnaseq/results/counts/Mov10_featurecounts.txt \
	  ~/ngs_course/rnaseq/results/STAR/bams/*bam
	  
> If you wanted to collect the information that is on the screen as the job runs, you can run it using the `2> filename` redirection, this type of redirection will collect all the information from the standard output into a file.

	**DO NOT RUN THIS**

	$ featurecounts -T 4 -s 2\ 
	  -a ~/ngs_course/unix_lesson/reference_data/chr1-hg19_genes.gtf \
	  -o ~/ngs_course/unix_lesson/rnaseq/results/counts/Mov10_featurecounts.txt \
	  ~/ngs_course/rnaseq/results/STAR/bams/*bam \
	  2> ~/ngs_course/unix_lesson/rnaseq/results/counts/Mov10_featurecounts.stdout

The output of this tool is 2 files, *a count matrix* and *a summary file* that tabulates how many the reads were "assigned"/counted and the reason they remained "unassigned". Lets take a look at the summary file:
	
	$ less results/counts/Mov10_featurecounts.txt.summary
	
Now lets look at the count matrix:
	
	$ less results/counts/Mov10_featurecounts.txt
	
There is information about the genomic coordinates and the length of the gene, we don't need this for the next step, so we are going to extract the columns that we are interested in.
	
	$ cut -f1,7,8,9,10,11,12 results/counts/Mov10_featurecounts.txt > results/counts/Mov10_featurecounts.Rmatrix.txt

The next step is to clean it up a little further by modifying the header line:
	
	$ vim results/counts/Mov10_featurecounts.Rmatrix.txt

### Note on counting PE data

For Paired-end data the bam file contains information about whether both read1 and read2 mapped and if they were at roughly the correct distance from each other, that is to say if they were "properly" paired. For most counting tools, only properly paired reads are considered by default, and each read pair is counted only once as a single "fragment". 

### Keeping track of read numbers

It is important to keep track of how many reads you started with and how many were counted as being associated with genes. This is important because it will help you pick out any obvious outlier, and it will also alert you early on to any issues with contamination and so on.

The things to keep track of for each sample are the following:
* number of raw reads
* number of reads left after trimming
* number of reads aligned to genome
* number of reads associated with genes 
 
We have an [Excel template](https://dl.dropboxusercontent.com/u/74036176/rna-seq_reads_template.xlsx) available, to get you started.

