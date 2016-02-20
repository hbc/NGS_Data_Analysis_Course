---
title: "Gene-level differential expression analysis using DESeq2"
author: "Meeta Mistry"
date: "Friday February 19, 2016"
---

Approximate time: 

## Learning Objectives 

* 


## DESeq2: Differential expression analysis

Now that the object is created we can run the `DESeq()` function. 

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

What just happened? Everything from normalization to linear modeling was carried out by the use of a single function, and the results of each step were inserted into the object that you initialized.

![deseq1](../img/deseq_obj2.png)



---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*