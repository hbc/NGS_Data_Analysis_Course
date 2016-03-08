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


module load seq/B

mkdir overlap_spp_macs2/

cat spp/Nanog_Rep1.narrowPeak spp/Nanog_Rep2.narrowPeak > overlap_spp_macs2/spp_Nanog.narrowPeak
sort -k1,1 -k2,2n overlap_spp_macs2/spp_Nanog.narrowPeak > overlap_spp_macs2/spp_Nanog_sorted.narrowPeak
bedtools merge -i overlap_spp_macs2/spp_Nanog_sorted.narrowPeak > overlap_spp_macs2/spp_Nanog_merged.bed 

cat macs2/Nanog-rep1_peaks.narrowPeak macs2/Nanog-rep2_peaks.narrowPeak > overlap_spp_macs2/macs2_Nanog.narrowPeak
sort -k1,1 -k2,2n overlap_spp_macs2/macs2_Nanog.narrowPeak > overlap_spp_macs2/macs2_Nanog_sorted.narrowPeak
bedtools merge -i overlap_spp_macs2/macs2_Nanog_sorted.narrowPeak > overlap_spp_macs2/macs2_Nanog_merged.bed 

bedtools intersect -a spp_Nanog_merged.bed -b macs2_Nanog_merged.bed -wo > Nanog_spp-macs_overlap.bed


module load seq/samtools

ll -htr *aln.bam

for fname in *aln.bam
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

