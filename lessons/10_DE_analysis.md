# Lesson

Differential expression analysis
===================

Learning Objectives:
-------------------
### What's the goal for this lesson?

* Use the count matrix as input to an R script for differential expression analysis
* Apply Unix commands to look at the results that are generated and extract relevant information
* Familiarize with various functional analysis tools



At the end of the workflow from the last lesson, our final end product was a count matrix. This is a matrix in which each row represents a gene (or feature) and each column corresponds to a sample. In our dataset, we have two sample classes (control and Mov10oe) and we want to assess the difference in expression between these groups on a gene-by-gene basis.

Intuitively, it would seem that since we know which samples belong to which group we could just compute a fold change for each gene and rank. Easy, right? Not exactly. The problem is, the gene expression that we are observing is not just a result of the differences between the groups that we are investigating, rather it is a measurement of the sum of many effects. In a given biological sample the transcriptional patterns will also be changing with respect to various extraneous factors; some that we are aware of (i.e demographic factors, batch information) and other noise that we cannot attribute to any particular source. So when we are looking for differentially expressed genes we need to try an control for these effects as best we can.


### Statistical models in R

Enter [R](https://www.r-project.org/), a software environment for statistical computing and graphics. R is widely used in the field of bioinformatics, amongst various other disciplines. It can be locally installed on almost all operating systems (and it's free!) with numerous packages available that help in increasing efficency of data handling, data manipulation and data analysis. Unfortunately, discussing the specifics about R is outside the scope of this course. However, we encourage you take a look at some of the R resources listed below if you are interested in learning more. 

R is a powerful language that can be very useful for NGS data analysis, and there are many popular tools for working with RNA-Seq count data. Some of the popular tools include [edgeR](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf), [DESeq2](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf), and [limma-voom](http://www.genomebiology.com/2014/15/2/R29). All of these tools use statistical modeling of the count data to test each gene against the null hypothesis and evaluate whether or not it is significantly differentially expressed. These methods account for within-group and between-group varaibility and are also flexible enough to allow for other coavriates you may want to account for. The details on how each tool works are described thoroughly within the vignettes.


### Setting up

We will be running an R script that uses the R package [edgeR](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) to identify differentially expressed genes. This package is available through the [Bioconductor](https://www.bioconductor.org/), a repository of packages for the analysis of high-throughput genomic data. 

While R is already installed on Orchestra for us to use, the package `edgeR` is not. To install the package we will need to first open up R. You can do this by simply typing `R` at the command prompt and pressing `Enter`. You are now in the R console:

![Rconsole](../img/R_screenshot.png)

To install the package `edgeR` you will need to copy/paste the lines of code below. The first line will source the Bioconductor installation script, and the second line uses the `biocLite()` function to install the package specified in quotations:
 	
	source("http://bioconductor.org/biocLite.R")
	biocLite("edgeR")

If this installed successfully, you can now exit R:

	q()


You should find yourself back at the shell command prompt. We will run the R script from the `results` directory, so let's navigate there and create a directory for the results of our analysis. We will call the directory `diffexpression`:

	cd results
	mkdir diffexpression

Now we need to copy over the required files. There are two files that we need, both are in different locations. First, let's copy over the script file:

	cp /groups/hbctraining/unix_oct2015_other/DE_script.R diffexpression/

The R script will require as input **1) your count matrix file** and **2) a metadata file**. The count matrix we generated in the last lesson and is in the `counts` directory. The metadata file is a tab-delimited file which contains any information associated with our samples. Each row corresponds to a sample and eaach column contains some information about each sample. Once you have the files copied, take a look at the metadata using `less`.

	cp counts/Mov10_rnaseq_counts_complete.txt diffexpression
	cp ~/unix_oct2015/other/Mov10_rnaseq_metadata.txt diffexpression


Ok, now we're all setup to run our R script! Let's run it from our `diffexpression` directory

	cd diffexpression
	Rscript DE_script.R Mov10_rnaseq_counts_complete.txt Mov10_rnaseq_metadata.txt 


### Resources for R
