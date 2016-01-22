---
title: "SAM/BAM file and assessing quality "
author: "Meeta Mistry, Bob Freeman"
date: "Wednesday, October 7, 2015"
---

Approximate time:

## Learning objectives



### Alignment file format: SAM/BAM

The output we requested from STAR is a BAM file, and by default returns a file in SAM format. BAM is a binary version of the SAM file, also known as Sequence Alignment Map format. The SAM file is a tab-delimited text file that contains information for each individual read and its alignment to the genome. The file begins with a header, which is optional, followed by an alignment section.  If present, the header must be prior to the alignments and starts with '@'. Each line that follows corresponds to alignment information for a read. Each alignment line has **11 mandatory fields for essential mapping information** and a variable number of fields for aligner specific information.


**Need a more in-depth description of the SAM file format**

These fields are described briefly below, but for more detailed information the paper by [Heng Li et al](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) is a good start.

![SAM](../img/SAM_file.png)


Let's take a quick look at our alignment. To do so we first convert our BAM file into SAM format using samtools and then pipe it to the `less` command. This allows us to look at the contents without having to write it to file (since we don't need a SAM file for downstream analyses). We first need to load the samtools module:

	module load seq/samtools/1.2

```
$ samtools view -h results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam | less

```
 
Scroll through the SAM file and see how the fields correspond to what we expected.

****

**Exercise:**
On examining the SAM file

***


---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
