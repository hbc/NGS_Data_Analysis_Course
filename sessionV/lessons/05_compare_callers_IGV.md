---
title: "Discover the overlap between peak callers and visualization with IGV"
author: "Radhika Khetani"
date: "Monday, March 7th, 2016"
---

Contributors: Radhika Khetani, 

Approximate time: 90 minutes

## Learning Objectives

* Learn how to perform coordinate-based analysis using bedtools
* Explore ChIP-Seq data in IGV using Encode-generated data for context


Bedtools description

Setting up

	$ bsub -Is -q interactive bash
	
	$ cd ~/ngs_course/chipseq/results/
	
	$ module load seq/BEDtools/2.23.0

	$ mkdir overlap_spp_macs2/
	
	$ cd overlap_spp_macs2/

Get overlap for Nanog
4 steps:
* combine replicates for each caller
* sort/re-order the combined file by coordinates for the combined files
* merge together any overlapping peaks (needs coordinate-sorted data) for each sorted file
* intersect the merged files and get overlapping peaks 

	$ cat ../spp/Nanog_Rep1.narrowPeak ../spp/Nanog_Rep2.narrowPeak > spp_Nanog.narrowPeak
	
	$ cat ../macs2/Nanog-rep1_peaks.narrowPeak ../macs2/Nanog-rep2_peaks.narrowPeak > macs2_Nanog.narrowPeak
	
	$ sort -k1,1 -k2,2n spp_Nanog.narrowPeak > spp_Nanog_sorted.narrowPeak
	
	$ sort -k1,1 -k2,2n macs2_Nanog.narrowPeak > macs2_Nanog_sorted.narrowPeak
	
	$ bedtools merge -h
	$ bedtools merge -i spp_Nanog_sorted.narrowPeak > spp_Nanog_merged.bed 

	$ bedtools merge -i macs2_Nanog_sorted.narrowPeak > macs2_Nanog_merged.bed 
	
	$ bedtools intersect -h
	$ bedtools intersect -a spp_Nanog_merged.bed -b macs2_Nanog_merged.bed -wo > Nanog_spp-macs_overlap.bed

	$ wc -l ../[sm]*/*Nanog*narrowPeak

	$ wc -l *Nanog*bed

Make .bai files for visualization

	cd ../bowtie2/
	
	$ module load seq/samtools

	$ ll -htr *aln.bam

	$ for fname in *aln.bam
		do
		samtools index $fname
		done

bring over with FileZilla:

bowtie2/H1hesc_Input_Rep1_chr12_aln.bam
bowtie2/H1hesc_Input_Rep1_chr12_aln.bam.bai
bowtie2/H1hesc_Nanog_Rep1_chr12_aln.bam
bowtie2/H1hesc_Nanog_Rep1_chr12_aln.bam.bai
spp/Nanog_Rep1.enrichment.estimates.wig
overlap_spp_macs2/Nanog_spp-macs_overlap.bed


IGV:
GRIN2B
SOX5
GPR19

investigate the blue mass on the graph

collapse the reads

bring in the nanog-rep1 peaks from encode

bring in broad peaks for H3K4me1 (mark for active/poised enhancers)

