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

We have already aligned the data, then sorted and cleaned it up for variant calling. In this module we will be calling variants using the program *[freebayes](https://github.com/ekg/freebayes#freebayes-a-haplotype-based-variant-detector)* ([Garrison and Marth, 2012)[http://arxiv.org/abs/1207.3907]), and filtering out low quality variants using *[vcftools](https://vcftools.github.io/index.html)*, a toolkit designed to work with Variant Call Format or VCF files.

<img src="../img/variant_calling_workflow_2.png" width="600">

## Freebayes

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




