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


## Variant Calling Workflow

The variant calling workflow begins with quality control and alignment, similar to the other NGS applications. Alignment is followed by alignment clean-up to prepare data for variant calling. Then, variant calling is performed, followed by filtering and annotation of the variant calls.

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




