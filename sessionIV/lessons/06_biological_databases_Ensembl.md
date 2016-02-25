---
title: "Biological Databases: NCBI"
author: "Mary Piper"
date: "Wednesday, October 7, 2015"
---

Contributors: Mary Piper

Approximate time: 3 hours

## Learning Objectives

* learn the pros/cons of different biological databases
* learn how to use features of biological databases to access information and data during an NGS analysis



## Intro to specific databases (Ensembl/Biomart) - 30 min. each

### The Ensembl Project
![ensembl_logo](../img/e_bang.png)

The [Ensembl](http://useast.ensembl.org/index.html) genome annotation system, developed jointly by the EBI and the Wellcome Trust Sanger Institute, has been used for the annotation, analysis and display of vertebrate genomes since 2000. All data are open source, i.e. freely available to the scientific community and is updated every 2-3 months.

Since 2009, the Ensembl site has been complemented by the creation of five new sites within [Ensembl genomes](http://ensemblgenomes.org), for [bacteria](http://ensemblgenomes.org/info/genomes?division=1), [protists](http://ensemblgenomes.org/info/genomes?division=5), [fungi](http://ensemblgenomes.org/info/genomes?division=3), [plants](http://ensemblgenomes.org/info/genomes?division=4) and [invertebrate metazoa](http://ensemblgenomes.org/info/genomes?division=2), enabling users to use a single collection of (interactive and programatic) interfaces for accessing and comparing genome-scale data from species of scientific interest from across the taxonomy. ***NOTE:** no annotations available through Ensembl for viral genomes.*

### Ensembl species and annotations
In the current release of Ensembl (83), over 80 vertebrate species are supported, with over half of these species being mammals. The dataset also includes the invertebrates yeast, *C. elegans*, and fruitfly to aid in more accurate generation of phylogenetic gene trees. *(Current species statistics for the  non-vertebrate databases are also [available](http://nar.oxfordjournals.org/content/44/D1/D574.full)).*

**All supported species have comprehensive, evidence-based gene annotations.** The "Gencode gene set" is used to create the Ensembl annotations and is made up of:

- Ensembl (automatically) annotated genes (using mRNA and protein sequences from UniProtKB and NCBI RefSeq) 
- Havana (manually) annotated genes
- Ensembl/Havana merges: transcripts that were identically annotated by both (reviewed annotations)

Gencode is the default gene set used by ENCODE, 1000 genomes and other major projects. 

![species_annotations](../img/species_annot.png)

A **selected set of genomes** includes additional data focused on *variation, comparative, evolutionary, functional and regulatory annotation*. The most advanced resources are provided for key species including **human, mouse, rat and zebrafish** reflecting the popularity and importance of these species in biomedical research.

![variation_species](../img/species_with_variation.png)


### Ensembl genome browser

#### Overview
Ensembl provides a genome browser that acts as a **single point of access to annotated genomes** for mainly vertebrate species. 

The browser can be used to easily access information at the genome, gene and protein level, such as gene sequence, splice variants, protein domains, genetic variation, homology, and regulatory elements. Ensembl imports genome sequences from consortia, which keeps the information consistent with many other bioinformatics projects. 

![ensembl_homepage](../img/ensembl_interface.png)

- **Searching Ensembl**:  Look for a gene, location, variant and more using the search box on the homepage or the box that is provided in the top right corner of any Ensembl page.

	- a gene name (for example, BRCA2) - best to use the official gene symbols ([HGNC](http://www.genenames.org))
	- a UniProt accession number (for example, P51587)
	- a disease name (for example, coronary heart disease)
	- a variation (for example, rs1223)
	- a location - a genomic region (for example, rat X:100000..200000)
	- a PDBe ID or a Gene Ontology (GO) term

	Most search results will take you to the appropriate Ensembl view through a results page. If you search using a location you will be directed straight to the location tab (this tab provides a view of a region of a genome).

- **Browse a Genome**: Choose your species of interest in this section. The drop down menu under 'All genomes' allows you to select from the full list. The *Ensembl Pre!* site contains new genomes (either new species to Ensembl, or updates in the reference assembly) that do not yet have an Ensembl gene set.  BLAST/BLAT is available for organisms in all Ensembl sites, including Pre!

- **Help**: There is a wealth of help and documentation in Ensembl if you are new to the browser. Video tutorials are provided, and printable pdfs with exercises. Custom data may be uploaded to Ensembl, or displayed directly by attaching a file by URL. 

- **News**: To find out what release you are working with, have a look at the news section of the homepage. If the current release is not the one you need, access archive sites to access previous versions, or releases, of Ensembl using the link on the lower right side.

#### Demo
Each species in Ensembl has its own home page, where you can find out who provided the genome sequence and which version of the genome assembly is represented.  

1. Click on `View full list of all Ensembl species` link. 
2. Click on the common name of your species of interest to go to the species homepage. We’ll click on `Human`.

	![ensembl_human](../img/ensembl_human.png)

	Within the human genome page, find the basic features:

	- Search bar for human information on gene, location, disease, etc.
	- News for the current human genome release
	- Information and sequences for the current human genome build
	- Links to example features in Ensembl
	- Guides on how to access information on comparative genomics, regulation, and variation
	
3. To find out more about genome assembly and the gene build, click on `More information and statistics`. Look over the information given in the statistics table.

	![stats](../img/ensembl_info.png)

4. Go back to the human genome page by clicking on the image. In the search bar type `mov10`.
5. From the search results select `MOV10 (Human Gene)`. The gene page for MOV10 should populate. 
	
	![gene_view](../img/ensembl_mov10_gene.png)
	
	The `Gene` page is organized as follows:
	
	- The side panel has **detailed gene information** displayed as a hierarchical tree by category.
	- The top of the page has a **gene overview**, giving a gene description, synonyms, location, and number of transcripts. 
	- This overview is followed by the **transcript table**. All transcripts based on any evidence are provided in the table. The transcripts are color-coded based on whether the transcript is protein-coding or non-coding, and the quality of evidence:
	
		- Gold: protein-coding transcripts that are Ensembl/Havana merges - essentially reviewed annotations with highest confidence
		- Red: protein-coding transcripts are less confidence
		- Blue: non-coding transcripts
		
		In addition to coloring, Ensembl also provides flags in the table for **"Transcript Support Levels"**, which highlight how well-supported or poorly-supported are the transcript models.
		
		Also provided in the table are the links to the **Consensus CoDing Sequence** sets (CCDS) for available transcripts. The CCDS is a consensus set of coding sequences established as a collaborative effort between NCBI, Ensembl, Vega, UniProt-SwissProt, and UCSC. 
		
		Note the Ensembl ID for MOV10 gene: ENSG00000155363. Ensembl uses the following format for naming:
	
		- ENSG###########	Ensembl Gene ID
		- ENST###########	Ensembl Transcript ID
		- ENSP###########	Ensembl Peptide ID
		- ENSE###########	Ensembl Exon ID
	
		For non-human species a suffix is added:
	
		- MUS (Mus musculus) for mouse ENSMUSG###
		- DAR (Danio rerio) for zebrafish: ENSDARG###
	

	- Below the transcript table is a summary section with links to external databases, followed by visualization of the transcripts.
	
	![transcript_vis](../img/ensembl_transcripts.png)
	
	**NOTE:**, in the visualization, the blue bar represents the genome contig, and transcripts above the bar are forward-stranded and those below are reverse.
	
6. Click on the MOV10-001 transcript ID, `ENST00000413052`. This should open a new tab entitled "Transcript: MOV10-001" with summary information for the transcript. In addition, the options on the side panel change.
		 
## Overview of features and interface without activities (can follow along if desired) - 10 min.
- basic or useful features of the database




### Use of biological databases in NGS analysis
Hands-on activities addressing important uses of biological databases during an NGS analysis

1. Hypothesis generation / exploration of genes of interest - Search for a gene and explore some of the info available, including basic info., sequence, isoform info., visualization, etc. - 10 min
2. Access to genome and gene annotation files - show where to find ftp sites - 5 min
3. Bringing data into an analysis - finding homologs, converting annotation ids, retreiving other NGS analysis data, etc. - 15-25 min.
4. Exploring NGS analysis results - population frequencies for variant calls, binding motifs, etc - only if have time


## Formatting	
The remaining sections are for formatting only. To maintain consistency throughout our lessons, we can format exercises and commands as follows.	

****
**Exercise**

1) Search for the sequence CTCAATGAAGAAATCTCTTAAAC in `Mov10_oe_1.subset.fq`.
In addition to finding the sequence, have your search also return
the name of the sequence.

2) Search for that sequence in all Mov10 replicate fastq files.
****

Let's try it out and put all the sequences that contain 'NNNNNNNNNN'
from all the files in to another file called `bad_reads.txt`.

`$ grep -B1 -A2 NNNNNNNNNN Mov10_oe_1.subset.fq > bad_reads.txt`
   
'-f1,3,4,5,7' means to cut these fields (columns) from the dataset.  

	chr1	exon	14362	14829	-
	chr1	exon	14970	15038	-
	chr1	exon	15796	15947	-
	chr1	exon	16607	16765	-
	chr1	exon	16858	17055	-

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
