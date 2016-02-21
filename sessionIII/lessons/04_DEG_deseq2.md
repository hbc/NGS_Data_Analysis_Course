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

### Results using a Wald test

To build a results table, we use the `results()` function on the `DESeqDataSet` object. By default, it will return to us the log2 fold changes and p-values for a Wald-test comparison of the last level over the first level. 

So in our case this would be control versus Mov1_0overexpression, and you can see that printed at the top of the output:

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

Conveniently, the default settings returned to us the Mov10_overexpression comparison to control, but it is usually best practice to specify which comparisons we are interested in looking at. This is especially useful when working with more than two factor levels. 

The comparisons are provided in the form of **contrasts**, in one of three different ways. In this lesson we will demonstrate the method that is most intuitive. By providing contrasts we are telling DESeq2 which coefficients to use for the hypothesis testing procedure. Let's take a look at the coefficients table (using `coef()`) to get an idea of what we have to choose from:

	View(coef(dds))

In total, we have four coefficients including the intercept. Each **coefficient** corresponds to a different factor level, **indicating the overall expression strength** of the gene. 

To specify the specific coeficients we are interested in, we need to provide the column names from the coefficents table as a list of 2 character vectors:

	contrast_kd <- list( "sampletypeMOV10_knockdown", "sampletypescontrol" )

**The order of the names, determines the direction of fold change that is reported.** The name provided in the second element is the level that is used to baseline. So for example, if we observe a fold change of -2 this would mean the gene expression is lower in Mov10_kd relative to the control. Pass the contrast as an argument to the `results()` function and take a look at the table that is returned:

	res_tableKD <- results(dds, contrast=contrast_kd)

***

**Exercise**

1. Create a contrasts vector for the Mov10_overexpression comparison to control.
2. Create a contrasts vector for the MOv10_overexpression comparison to *all other samples*.

*** 



### Exploring the results table

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*