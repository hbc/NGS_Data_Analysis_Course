---
title: "Counting reads"
author: "Meeta Mistry, Bob Freeman"
date: "Wednesday, October 7, 2015"
---

Approximate time: 

## Learning Objectives:

*
* 


### Counting reads
Once we have our reads aligned to the genome, the next step is to count how many reads have been mapped to each gene. Counting is done with a tool called [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html). The input files required for counting include the BAM file and an associated gene annotation file in GTF format. htseq-count works by taking the alignment coordinates for each read and cross-referencing that to the coordinates for features described in the GTF. Most commonly a feature is considered to be a gene, which is the union of all exons (which is a feature type) that map to that gene. There is no minimum overlap to determine whether or not a read is counted for a particular gene, rather it is the mode that the user chooses. 


We will be using the 'union' mode as it is default and most commonly used. To find out more on the different modes and how they affect your output, take a look at the [manual](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)
	
	$ featureCounts 

	  Version 1.4.4

	  Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 



Let's start by creating a directory for the output:

	$ mkdir results/counts

	$ export PATH=/opt/bcbio/local/bin:$PATH


In it's most basic form the htseq command requires only the BAM file and the GTF file. We will add in a few additional parameters including `--format` to indicate BAM file, and `--stranded reverse` to specify that we have a stranded library created via the dUTP method. By default htseq-count will _ignore any reads that map to multiple locations_ on the genome. This results in undercounting but also helps reduce false positives. While multi-mappers are a feature that cannot be modified, there is a parameter that allows the user to filter reads by specifying a minimum alignment quality. 

You will notice at the end of the command we have added a redirection symbol. Since htseq-count outputs results to screen, we need to re-direct it to file.

	featureCounts -T 6 -a ~/ngs_course/unix_lesson/reference_data/chr1-hg19_genes.gtf \
		-o ~/ngs_course/unix_lesson/rnaseq/results/counts/counts.txt \
		~/ngs_course/unix_lesson/rnaseq/results/STAR/*bam


#### Exercise
Take a look at the end of the file using the `tail` command. You should see a summary of how the reads were classified. 

1. How many reads were assigned as no_feature? Why would they be classified this way?
2. How many reads were found to map to multiple locations?






