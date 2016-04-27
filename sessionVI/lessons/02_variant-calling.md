---
title: "Variant calling with Freebayes"
author: "Radhika Khetani, Meeta Mistry"
date: "Thursday, April 28th, 2016"
---

Approximate time: 45 minutes

## Learning Objectives:

* Call variants with Freebayes
* Get familiar with the Variant Call Format (VCF)
* Use vcftools to perform some simple filtering on the variants in the VCF file


## Variant Calling

We have already aligned the data, then sorted and cleaned it up for variant calling. 

<img src="../img/variant_calling_workflow_2.png" width="600">

Some of the more popular tools for calling variants include [SAMtools](http://samtools.sourceforge.net/mpileup.shtml), [the GATK suite](https://www.broadinstitute.org/gatk/about/) and [FreeBayes](https://github.com/ekg/freebayes#freebayes-a-haplotype-based-variant-detector) ([Garrison and Marth, 2012](http://arxiv.org/abs/1207.3907)). While it can be useful to work through the [GATK Best Practices](https://www.broadinstitute.org/gatk/guide/best-practices.php) we will be using FreeBayes in this module as it is just as sensitive and precise, but has no license restrictions. After calling variants, we will filter out low quality variants using *[vcftools](https://vcftools.github.io/index.html)*, a toolkit designed to work with Variant Call Format or VCF files.

## Freebayes

*FreeBayes* is a **haplotype-based** variant detector and is a great tool for calling variants from a population. To understand how Freebayes works, let's take a look at their [user manual and guide](https://github.com/ekg/freebayes#freebayes-a-haplotype-based-variant-detector).

<img src="../img/freebayes_1.png" width="600">

	$ mkdir variants
	$ cd variants/
	$ which freebayes
	$ freebayes -f ../../data/reference_data/chr20.fa \
	../bwa/na12878_sorted_marked.bam > na12878.vcf


## Variant Call Format (VCF)

	$ less na12878.vcf

## Filtering VCFs

	$ module avail seq/vcf
	$ module load seq/vcftools/0.1.12
	$ vcftools --vcf na12878.vcf --minQ 20 --recode --recode-INFO-all \
	--out na12878_q20  
	
***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*




