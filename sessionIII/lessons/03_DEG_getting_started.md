---
title: "Gene-level differential expression analysis using DESeq2"
author: "Meeta Mistry"
date: "Friday February 19, 2016"
---

Approximate time: 

## Learning Objectives 

* 


## Differential expression analysis
Thus far, we have described different strategies for RNA-Seq data (i.e. de novo transcriptome assembly, transcript discovery) but arguably the most common use for transcriptome data is to to search for differentially expressed genes. Finding genes that are differentially expressed between conditions is an integral part of understanding the molecular basis of phenotypic variation [[Soneson and Dleorenzi, 2013](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-91)].

There are a number of software packages that have been developed for differential expression analysis of RNA-seq data, and new methods are continuously being presented. Many studies describing comparisons between these methods show that while there is some agreement, there is also much variability. **Additionally, there is no one method that performs optimally under all conditions.**


![deg1](../img/deg_methods1.png) 

![deg1](../img/deg_methods2.png) 


In the next few lessons, we will walk you through an end-to-end gene-level RNA-seq differential expression workflow using various R packages. We will start with a count matrix and perform exploratory data analysis for quality assessment and to explore the relationship between samples, perform differential expression analysis, and visually explore the results.

## Setting up

Let's get started by opening up RStudio and setting up a new project for this analysis. 

1. Go to the `File` menu and select `New Project`.
2. In the `New Project` window, choose `New Directory`. Then, choose `Empty Project`. Name your new directory `DEanalysis` and then "Create the project as subdirectory of:" the Desktop (or location of your choice).
3. After your project is completed, it should automatically open in RStudio. 

To check whether or not you are in the correct working directory, use `getwd()`. The path `Destktop/DEanalysis` shoudl be returned to you in the console. Within your working directory use the `New folder` button in the bottom righ panel to create three new directories: `data`, `meta` and `results`. Remember the key to a good analysis is keeping organized from the start!

Go to the `File` menu and select `New File`, and select `R Script`. This should open up a script editor in the top left hand corner. This is where we will be typing and saving all commands required for this analysis. In the script editor type in a header line:

```
## Gene-level differential expression analysis using DESeq2
## NGS Data Analysis 2016
```

Now save the file as `de_script.R`. When finished your working directory should now look something like this:

![setup](../img/settingup.png)

Finally, we need to grab the files that we will be working with for the analysis. For each of the files below we have provided a link. Right click on the link, and "Save link as ...". The counts matrix will need to be put inside the `data` directory, and the metadata table in your `meta` directory.

* [Full counts matrix](https://raw.githubusercontent.com/hbc/NGS_Data_Analysis_Course/master/sessionIII/data/Mov10_full_counts.txt)
* [Full metadata table](https://raw.githubusercontent.com/hbc/NGS_Data_Analysis_Course/master/sessionIII/data/Mov10_full_meta.txt)

## Loading libraries

For this analysis we will be using several R packages, some which have been installed from CRAN and others from Bioconductor. To use these packages (and the functions contained within them), we need to **load the libraries.** Add the following to your script and don't forget to comment liberally!

```
## Setup
### Bioconductor and CRAN libraries used
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
```

## Loading data

To load the data into our current environment, we will be using the `read.table` function. We need to provide the path to each file and also specify arguments to let R know that we have a header (`header = T`) and the first column is our row names (`row.names =1`). By default the function expects tab-delimited files, which is what we have.

```
## Load in data
data <- read.table(file.path(dataDir, 'Mov10_full_counts.txt'), header=T, 
row.names=1, as.is=T) 
meta <- read.table(file.path(metaDir, 'Mov10_full_meta.txt'), header=T, 
row.names=1)
```

Use `class()` to inspect our data and make sure we are working with data frames:

```
### Check classes of the data we just brought in
class(data)
class(meta)
```

As a sanity check we should also make sure that we have sample names that match between the two files, and that the samples are in the right order.

	### Check that sample names match in both files
	all(names(data) %in% rownames(meta))
	all(names(data) == rownames(meta))
	
***

**Exercise**	

1. Suppose we had sample names matching in the counta matrix and metadata file, but they were out of order. Write the line(s) of code required to create a new matrix with columns ordered such that they were identical to the row names of the metadata.

*** 


## `DESeq2DataSet`

Now that we have all the required files and libraries loaded we are ready to begin with the exploratory part of our analysis. 

Amongst the methods that work on count data directly, a popular tool is [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). The methodology of DESeq2 (as described in the lecture), has high sensitivity and precision, while controlling the false positive rate. 

The first thing we need to do is create `DESeqDataSet` object. Bioconductor software packages often define and use a custom class for storing data that makes sure that all the needed 'data slots' are consistently provided and fulfill the requirements. These objects are similar to `lists` in that the `data slots` are analagous to components as they store a number of different types of data structures. These objects are **different from lists** in that the slots are designated for specific information and access to that information (i.e. selecting data from the object) is done using object-specific functions.

Let's start by creating the `DESeqDataSet` object and then we can talk a bit more about what is stored inside it. To create the object we will need the **count matrix** and the **metadata** table as input. We will also need to specify a **design formula**. The design formula specifies the column(s) in the metadata table and hwo they should be used in the analysis. For our data set we only have one column we are interested in, that is `~sampletype`. This column has three factor levels, which tells DESeq2 that for each gene we want to evaluate gene expression change with respect to these different levels.


	## Create DESeq2Dataset object
	dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
					design = ~ sampletype)


![deseq1](../img/deseq_obj1.png)


You can use DESeq-specific functions to access the different slots and retrieve information if you wish. For example, suppose we wanted the original count matrix we would use `counts` (*Note: we nested it within the `View` function so that rather than getting printed in the console we can see it in the script editor*) :

	View(counts())

As we go through the workflow we will use the relevant functions to check what information gets stored inside our object.

## Quality assessment and exploratory analysis



	### Transform counts for data visualization
	rld <- rlog(dds, blind=TRUE)
	rld_mat <- assay(rld) 


---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*