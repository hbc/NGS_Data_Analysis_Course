---
title: "Biological Databases"
author: "Mary Piper"
date: "Wednesday, October 7, 2015"
---

Contributors: Mary Piper

Approximate time:  minutes

## Learning Objectives

* learn the pros/cons of different biological databases
* learn how to use biological databases during an NGS analysis

## Intro to biological databases - 15-20 min.

This section should be covered in a lecture-based format. It will include:

- primary versus derivative databases
- discussion of features of primary databases, such as how often updated, etc.
- no hands-on exercises
- no discussion of Ensembl/NCBI/UCSC interfaces or features


## Intro to specific databases (Ensembl/NCBI/UCSC) - 50 min. each

### Overview of features and interface without activities (can follow along if desired) - 10 min.
- basic or useful features of the database

### Use of biological databases in NGS analysis
Hands-on activities addressing important uses of biological databases during an NGS analysis

1. Hypothesis generation / exploration of genes of interest - Search for a gene and explore some of the info available, including basic info., sequence, isoform info., visualization, etc. - 10 min
2. Access to genome and gene annotation files - show where to find ftp sites - 5 min
3. Bringing data into an analysis - finding homologs, converting annotation ids, retreiving other NGS analysis data, etc. - 15-25 min.
4. Exploring NGS analysis results - population frequencies for variant calls, binding motifs, etc - only if have time


	
	

****
**Exercise**

1) Search for the sequence CTCAATGAAGAAATCTCTTAAAC in `Mov10_oe_1.subset.fq`.
In addition to finding the sequence, have your search also return
the name of the sequence.

2) Search for that sequence in all Mov10 replicate fastq files.
****

## Redirection

Let's try it out and put all the sequences that contain 'NNNNNNNNNN'
from all the files in to another file called `bad_reads.txt`.

`$ grep -B1 -A2 NNNNNNNNNN Mov10_oe_1.subset.fq > bad_reads.txt`
   
'-f1,3,4,5,7' means to cut these fields (columns) from the dataset.  

	chr1	exon	14362	14829	-
	chr1	exon	14970	15038	-
	chr1	exon	15796	15947	-
	chr1	exon	16607	16765	-
	chr1	exon	16858	17055	-


