---
title: "QC + Alignment with BWA"
author: "Meeta Mistry, Mary Piper"
date: "Wednesday, March 23, 2016"
---

Approximate time: 90 minutes

## Learning Objectives:

* Exploration of the variant calling workflow
* Understanding the alignment method BWA utilizes to align sequence reads to the reference genome
* Choosing appropriate BWA alignment parameters for our dataset

## Variant Calling Workflow

![var_calling_workflow](../img/variant_calling_workflow.png)

The variant calling workflow begins with quality control and alignment, similar to the other NGS applications. Alignment is followed by alignment cleanup to prepare data for variant calling. Then, variant calling is performed, followed by variant call filtering and annotation of the variant calls.

## Read Alignment

Choice of alignment tool is often determined by the type of NGS application being conducted. We have previously used STAR for RNA-Seq data because it is fast and optimized for aligning spliced reads. For ChIP-Seq we used Bowtie2 to align the reads because it is fast and accurate. For variant calling we will use [BWA (Burrows-Wheeler Aligner)](http://bio-bwa.sourceforge.net) for alignment. 

BWA is generally slower than Bowtie2 with similar sensitivity and both tools can perform gapped alignment for the identification of indels and can effectively map paired-end reads. However, BWA is a bit more accurate and provides information on which alignments are trustworthy. Small numbers of bad alignments can result in many false variant calls, so accuracy is paramount when performing variant calling.

### BWA modes

Depending on read length, BWA has different modes optimized for different sequence lengths:

- **BWA-backtrack:** designed for Illumina sequence reads up to 100bp (3-step)

- **BWA-SW:** longer sequences ranged from 70bp to 1Mbp, long-read support and split alignment

- **BWA-MEM:** share similar features to BWA-SW, but BWA-MEM is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

### Aligning reads with BWA

#### Set-up


For all the algorithms, BWA first needs to construct the FM-index for the reference genome  .




