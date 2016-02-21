---
title: "Gene-level differential expression analysis using DESeq2"
author: "Meeta Mistry"
date: "Friday February 19, 2016"
---

Approximate time: 

## Learning Objectives 

* 


## DESeq2: Differential expression analysis

### Getting setup

Let's get started by opening RStudio and opening up the project that we created yesterday. 

1. Go to the File menu and select 'Open project ...'
2. Navigate to `~/Desktop/DEanalysis/` and double click on the `DEanalysis.Rproj` file

You should see your environment become populated with all of the variables created last lesson. The only thing that we will need to do is reload the required libraries:

```
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
```


### Running DESeq2

To run the differential expression pipeline on the raw counts in DESeq2, we use a **single call to the function `DESeq()`**. The required input is the `DESeqDataSet` object that we created in the last lesson. By re-assigning the results of the function back to the same variable name, we can continue to fill in the `slots` in our `DESeqDataSet` object.

	##Run analysis
	dds <- DESeq(dds)
 
This function will print out a message for the various steps it performs. 

```
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
``` 

What just happened? Remember our workflow:

![workflow](../img/deseq_workflow1.png)


Everything from normalization to linear modeling was carried out by the use of a single function! The results of each step were inserted into the object that you initialized.

![deseq1](../img/deseq_obj2.png)

But wait, you may be asking yourself shouldn't the QC assessment have been performed post-normalization? The answer is, yes it was.

When you ran the `rlog()` function the normalization was performed. And so if we compare the `sizeFactors` value in both objects you will find that the values are identical.

	sizeFactors(dds)
	rld$sizeFactors

> *NOTE:* The two objects require different ways in which we can access the information stored inside. When in doubt, use `class()` to find out what type of data structure you are working with. Knowing this information is key to finding ways of extracting infomration from the object.

### Extracting results using a Wald test

To build a results table, we use the `results()` function on the `DESeqDataSet` object. By default, it will return to us the log2 fold changes and p-values for a Wald-test comparison of the last level over the first level. 

So in our case this would be control versus Mov1_0overexpression, and you can see that printed at the top of the output:

	## Extract results
	res_tableOE <- results(dds)
	head(res_tableOE)

```
log2 fold change (MAP): sampletype MOV10_overexpression vs control 
Wald test p-value: sampletype MOV10_overexpression vs control 
DataFrame with 6 rows and 6 columns
               baseMean log2FoldChange      lfcSE       stat    pvalue       padj
              <numeric>      <numeric>  <numeric>  <numeric> <numeric>  <numeric>
1/2-SBSRNA4  45.6520399     0.26976764 0.18775752  1.4367874 0.1507784 0.25242910
A1BG         61.0931017     0.20999700 0.17315013  1.2128030 0.2252051 0.34444163
A1BG-AS1    175.6658069    -0.05197768 0.12366259 -0.4203185 0.6742528 0.77216278
A1CF          0.2376919     0.02237286 0.04577046  0.4888056 0.6249793         NA
A2LD1        89.6179845     0.34598540 0.15901426  2.1758136 0.0295692 0.06725157
A2M           5.8600841    -0.27850841 0.18051805 -1.5428286 0.1228724 0.21489067
```

Let's go through some of the columns in the results table to get a better idea of what we are looking at. To extract information regarding the meaning of each column we can use `mcols()`:

	mcols(results_wald, use.names=T)

* `baseMean`: mean of normalized counts for all samples
* `log2FoldChange`: log2 fold change
* `lfcSE`: standard error
* `stat`: Wald statistic
* `pvalue`: Wald test p-value
* `padj`: BH adjusted p-values
 


#### Contrasts

Conveniently, the default settings returned to us the Mov10_overexpression comparison to control, but it is usually best practice to specify **which comparisons we are interested in** looking at. This is especially useful when working with more than two factor levels or more complex designs.

The comparisons are provided in the form of **contrasts**, in one of three different ways. In this lesson we will demonstrate the method that is most intuitive. By providing contrasts we are telling DESeq2 which coefficients to use for the hypothesis testing procedure. Let's take a look at the coefficients table (using `coef()`) to get an idea of what we have to choose from:

	View(coef(dds))

In total, we have four coefficients including the intercept. Each **coefficient** corresponds to a different factor level, **indicating the overall expression strength** of the gene. 

To specify the specific coeficients we are interested in, we need to provide the column names from the coefficents table as a list of 2 character vectors:

	## Define contrasts
	contrast_kd <- list( "sampletypeMOV10_knockdown", "sampletypescontrol" )

**The order of the names, determines the direction of fold change that is reported.** The name provided in the second element is the level that is used to baseline. So for example, if we observe a fold change of -2 this would mean the gene expression is lower in Mov10_kd relative to the control. Pass the contrast as an argument to the `results()` function and take a look at the table that is returned:

	res_tableKD <- results(dds, contrast=contrast_kd)

***

**Exercise**

1. Create a contrasts vector for the Mov10_overexpression comparison to control.
2. Create a contrasts vector for the Mov10_overexpression comparison to *all other samples*.

*** 



### Summarizing the results table

To summarize the results table, a handy function in DESeq2 is `summary()`. Confusingly it has the same name as the function used to inspect data frames. This function when called with a DESeq results table as input, will summarize the results at a given FDR threshold. 

	## Suammrize results
	summary(res_tableOE)
	

```  
out of 19748 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 3657, 19% 
LFC < 0 (down)   : 3897, 20% 
outliers [1]     : 0, 0% 
low counts [2]   : 3912, 20% 
(mean count < 4)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

In addition to the number of genes up- and down-regulated at and FDR < 0.1, the function also reports the number of genes that were tested (genes with non-zero total read count), and the number of genes not included in multiple test correction due to a low mean count (which in our case is < 4).

The default FDR threshold is set to `alpha = 0.1`, which is quite liberal. Let's try changing that to `0.05` -- *how many genes are we left with*?

The FDR threshold on it's own doesn't appear to be reducing the number of significant genes. With large significant gene lists it can be hard to extract meaningful biological relevance. To help increase stringency, one can also add a fold change threshold. The `summary()` function doesn't have an argument for fold change threshold, but instead we can use the base R function `subset()`.

Let's first create variables that contain our threshold criteria:

	### Set thresholds
	padj.cutoff <- 0.05
	lfc.cutoff <- 1

The `lfc.cutoff` is set to 1; remember that we are working with log2 fold changes so this translates to an actual fold change of 2 which is pretty reasonable. Now let's setup our **`subset()` function nested within the `summary()` function**. Start building from the inside out:

	subset(res_tableOE)

We need to add our selection criteria. The first is our FDR threshold:

	subset(res_tableOE, padj < padj.cutoff)

Now let's add in the log2 fold change criteria. Because we want both up- and down-regulated genes we will use the absolute value of the fold change using the `abs(log2FoldChange)` function:

	subset(res_tableOE, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

Now, finally we will put all of that inside the `summary()` function. This is a fast way of getting overall statistics and deciding whether our threshold is still too liberal or perhaps overly stringent.

	summary(subset(res_tableOE, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff))


Does this reduce our results? How many genes are up-regulated and down-regulated at this new threshold?

We should have a total of 86 genes that are significantly differentially expressed. To denote these genes as significant we can add a column in our results table. The vector will be a logical vector, where `TRUE` means the gene passes our threshold and `FALSE` means it fails.

	res_tableOE$threshold <- as.logical(res_tableOE$padj < padj.cutoff & 
                   abs(res_tableOE$log2FoldChange) > lfc.cutoff)

Now we can easily check how many genes are significant by using the `which()` function:

	length(which(res_tableOE$threshold))

***

**Exercise**

1. Explore the results table for the Mov10_knockdown comparison to control. How many genes are differntially expressed using the default thresholds?
2. Using the same thresholds set above (`padj.cutoff <- 0.05` and `lfc.cutoff <- 1`), report he number of genes that are up- and down-regulated in Mov10_knockdown compared to control.
3. Add a new column called `threshold` to the `res_tableKD` which contains a logical vector denoting genes as being differentially expressed or not

*** 


## Visualizing the results





## Exporting significant gene lists

The next step in our workflow is interpretation of gene lists using various tools for functional analysis. Depending on the tool you choose to use downstream, you will require different information from the results table as input. To be safe it is wise to keep atleast one copy of the full results table with relevant information. 


First, let's sort the results file by adjusted p-value:
	
	### Sort the results tables
	res_tableOE_sorted <- res_tableOE[order(res_tableOE$padj), ]
	res_tableKD_sorted <- res_tableKD[order(res_tableKD$padj), ]
	
Now we can use the `write.table()` function to write them to file:

	### Write sorted results to file
	write.table(res_tableOE_sorted, file="results/results_OE_sortedPval.txt", sep="\t", quote=F, col.names=NA)
	
	write.table(res_tableKD_sorted, file="results/results_KD_sortedPval.txt", sep="\t", quote=F, col.names=NA)

One of the tools we will be using for functional analysis will require only the gene names of the significant genes, but ordered by adjusted p-value. We can easily create this here since our results files are already sorted.

	### Get significant genes
	sigOE <- row.names(res_tableOE_sorted)[which(res_tableOE_sorted$threshold)]
	sigKD <- row.names(res_tableKD_sorted)[which(res_tableKD_sorted$threshold)]
	
To write these lists to file we will use the `write()` function which will write the contents to file on single line, or if `ncol` is specified, into a certain number of columns:

	### Write genes to file
	write(sigOE, file="results/Mov10_oe_logFC_1_pVal_0.05.txt", ncol=1)
	write(sigKD, file="results/Mov10_kd_logFC_1_pVal_0.05.txt", ncol=1)
	
## Saving the project

Now we are set up for functional analysis of our gene lists. Make sure you save your R session as you quit RStudio to your DEanalysis project, so you don't lose all your work from this DE analysis module!

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*