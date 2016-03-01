---
title: "Biological Databases: NCBI"
author: "Mary Piper"
date: "Wednesday, October 7, 2015"
---

Contributors: Meeta Mistry

Approximate time: 30 minutes

## Learning Objectives

* learn how to use features of the NCBI database and genome browser to access information and data during an NGS analysis

## Intro to NCBI

The [National Center for Biotechnology Information (NCBI)](http://www.ncbi.nlm.nih.gov/) is part of the United States National Library of Medicine (NLM), a branch of the National Institutes of Health. The NCBI is located in Bethesda, Maryland and was founded in 1988 through legislation sponsored by Senator Claude Pepper.

<img src="../img/ncbi_screenshot.png" width="800">

The website acts as a portal to a vast amount of biological information. There a number of databases and analysis available at your fingertips which allow you to explore your biological question. You can browse data, submit, download and even analyze data.

NCBI is most popular for its PubMed resource for searching the biological literature and BLAST for sequence similarity searching (shown in the 'Popular Resources' panel on the right). However, the functionality of NCBI goes far beyond those tools. 


### Searching NCBI

<img src="../img/ncbi_screenshot3.png" width="400">

There are many ways of performing a search in NCBI. You can simply type in a gene/gene product or keyword into the search box and search all databases for a hit, and then drill down to a more specific result from there. Alternatively, you can narrow down your search from the start the start.


### Building a search using filters

Let's start with just a simple gene search using the search box and selecting all databases.

> Type in Mov10 to the search box and press 'Enter`.

You should find there are may hits returned. MOv10 appears in 27 different databases. We can narrow this search down if we were interested only in gene information.

> Click on Gene within the Genes section of the results. How many entries are we now reduced to?

There are still ~1000 genes returned. Take a look at the hint box at the top of our results:

<img src="../img/hint-box1.png" width="600">

Since this applies to us, we can click on the link and reduce our results further.

> Take a look at the left-hand panel and other filters that are available to you. Activate one or two and see how it changes your results. 

On the right hand-side you should see how your search was "built" in the 'search details' box.

	Mov10[sym] AND (gene_nucleotide_pos[filter] AND alive[prop])

Rather than using the graphical user interface to filter, we can also build a search from scratch:

<img src="../img/search_builder.png" width="600">

> Try building the same search using filters, or select it from your history.

Take a look at the gene record for Mov10. There is a Table of Contents  on the left to help navigate through the wealth of information.

<img src="../img/gene-card.png", width="500">

There is also Mapview of the gene incorporated into the card, this is NCBI's genom browser. It is very complicated to interpret, but if you want to take a stab at it, this [graphical view legend](http://www.ncbi.nlm.nih.gov/tools/sviewer/legends/) is very helpful. 


### Downloading data from NCBI

There are many ways to download data, depending on where you are in NCBI. On the main page, there is an icon which will direct you to the FTP site, and you can search the folders to find what is appropriate.

<img src="../img/ncbi_download.png" width="200">

Alternatively, on the left-hand panel of the main screen, there are resources listed by category. This allows you easily navigate to what is relevant. 

1. Click on the "Resources List (A-Z)", and you should see resources listed alphabetically.

<img src="../img/az-ncbi.png" width="400">

2. Scroll down to the "F" section and you will see a list of links to FTP sites for each of the different biological data types that are available for download. 

<img src="../img/ftp_ncbi.png" width="200">

The most common use for genomic data is the genome data, as you have already encountered in this course. But other data types can also be useful for studies that involve data integration. 


***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
