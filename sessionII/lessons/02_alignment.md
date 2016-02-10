---
title: "Alignment with STAR"
author: "Meeta Mistry, Bob Freeman, Mary Piper"
date: "Wednesday, February 10, 2016"
---

Approximate time:

## Learning Objectives:

* Understanding the alignment method STAR utilizes to align sequence reads to the reference genome
* Identifying the intricacies of alignment tools used in NGS analysis (parameters, usage, etc)
* Choosing appropriate STAR alignment parameters for our dataset
* Running STAR on multiple samples

## Read Alignment

The alignment process consists of choosing an appropriate reference genome to map our reads against and performing the read alignment using one of several splice-aware alignment tools such as [STAR](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635) or [TopHat2](https://ccb.jhu.edu/software/tophat/index.shtml). The choice of aligner is a personal preference and also dependent on the computational resources that are available to you.

## STAR Aligner

To determine where on the human genome our reads originated from, we will align our reads to the reference genome using STAR (Spliced Transcripts Alignment to a Reference). STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping, and utilizes a novel strategy for spliced alignments. 

### STAR Alignment Strategy

STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed (but also requires quite a bit of memory). The algorithm achieves this highly efficient mapping by performing a two-step process:

1. Seed search
2. 
	

### Running STAR

#### Set-up

To get started with this lesson, start an interactive session with 6 cores:

```
$ bsub -Is -n 6 -q interactive bash	
```

You should have a directory tree setup similar to that shown below. it is best practice to have all files you intend on using for your workflow present within the same directory. In our case, we have our original FASTQ files and post-trimming data generated in the previous section. We also have all reference data files that will be used in downstream analyses.

```
rnaseq
	├── data
	│   ├── reference_data
	│   │   └── chr1.fa
	│   │   └── chr1-hg19_genes.gtf
 	|   ├── untrimmed_fastq
	│   │   
	│   └── trimmed_fastq
	│       ├── Irrel_kd_1.qualtrim25.minlen35.fq
	│       ├── Irrel_kd_2.qualtrim25.minlen35.fq
	│       ├── Irrel_kd_3.qualtrim25.minlen35.fq
	│       ├── Mov10_oe_1.qualtrim25.minlen35.fq
	│       ├── Mov10_oe_2.qualtrim25.minlen35.fq
	│       └── Mov10_oe_3.qualtrim25.minlen35.fq
	|
	├── meta
	├── results
	└── docs
```


Change into the rnaseq folder. This will be our working directory:

```
$ cd ngs_course/rnaseq

```

Let's also load up the STAR module: 

```
$ module load seq/STAR/2.4.0j
```

Create an output directory for our alignment files:

```bash
$ mkdir results/STAR

```

For now, we're going to work on just one sample to set up our workflow. To start we will use the trimmed first replicate in the Mov10 over-expression group, `Mov10_oe_1.subset.fq.qualtrim25.minlen35.fq`.


**NOTE: if you did not follow the last section, please execute the following command:** (this will copy over the required files into your home directory.)

```bash

$ cp -r /groups/hbctraining/unix_oct2015_other/trimmed_fastq data/ .

```

#### Aligning Reads


Aligning reads using STAR is a two step process:   

1. Create a genome index 
2. Map reads to the genome

> A quick note on shared databases for human and other commonly used model organisms. The Orchestra cluster has a designated directory at `/groups/shared_databases/` in which there are files that can be accessed by any user. These files contain, but are not limited to, genome indices for various tools, reference sequences, tool specific data, and data from public databasese such as NCBI and PDB. So when using a tool and requires a reference of sorts, it is worth taking a quick look here because chances are it's already been taken care of for you. 

##### Creating a genome index

Indexing of the reference genome has already been done for you. **You do not need to run this code**. For this step you need to provide a reference genome and an annotation file. For this workshop we are using reads that originate from a small subsection of chromosome 1 (~300,000 reads) and so we are using only chr1 as the reference genome, and have provided the appropriate indices. Depending on the size of your genome, this can take a while. 

The basic options to **generate genome indices** using STAR as follows:


* `--runThreadN`: number of threads
* `--runMode`: genomeGenerate mode
* `--genomeDir`: /path/to/store/genome_indices
* `--genomeFastaFiles`: /path/to/FASTA_file 
* `--sjdbGTFfile`: /path/to/GTF_file
* `--sjdbOverhang`: readlength -1

```
** Do not run this**
STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles chr1.fa --sjdbGTFfile chr1-hg19_genes.gtf --sjdbOverhang 99

```

##### Aligning reads

The basic options for aligning reads to the genome using STAR is as follows:

* `--runThreadN`: number of threads
* `--readFilesIn`: /path/to/FASTQ_file
* `--genomeDir`: /path/to/genome_indices
* `--outFileNamePrefix`: prefix for all output files


More details on STAR and its functionality can be found in the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), we encourage you to peruse through to get familiar with all available options.

We can access the software by simply using the STAR command followed by the basic parameters described above and any additional parameters. The full command is provided below for you to copy paste into your terminal. If you want to manually enter the command, it is advisable to first type out the full command in a text editor (i.e. [Sublime Text](http://www.sublimetext.com/) or [Notepad++](https://notepad-plus-plus.org/)) on your local machine and then copy paste into the terminal. This will make it easier to catch typos and make appropriate changes. 

Below, we first describe some the extra parameters we have added.

Advanced parameters:

* `--outFilterMultimapNmax`: max number of multiple alignments allowed for a read
* `--outSAMstrandField`: compatability with Cufflinks (for transcriptome assembly)
* `--outReadsUnmapped`: file format for unmapped reads
* `--outSAMtype`: output filetype (SAM default)
* `--outSAMUnmapped`: what to do with unmapped reads
* `--outSAMattributes`: specify SAM attributes in output file


```
STAR --runThreadN 6 --genomeDir /groups/hbctraining/unix_oct2015_other/reference_STAR \
--readFilesIn data/trimmed_fastq/Mov10_oe_1.qualtrim25.minlen35.fq \ 
--outFileNamePrefix results/STAR/Mov10_oe_1_ \
--outFilterMultimapNmax 10 \
--outSAMstrandField intronMotif \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI NM MD AS
```

****
**Exercise**

How many files do you see in your output directory? Using the `less` command take a look at `Mov10_oe_1_Log.final.out` and answer the following questions:  

1. How many reads are uniquely mapped?
2. How many reads map to more than 10 locations on the genome?
3. How many reads are unmapped due to read length?   

***



### Running STAR on multiple samples

Create a job submission script to run STAR on all FASTQ files

```
add script here 
with LSF directives

```

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). All instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*









