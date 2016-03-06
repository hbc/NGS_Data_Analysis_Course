---
title: "Differenial Peak calling MACS2"
author: "Meeta Mistry"
date: "Sunday, March 6th, 2016"
---

Contributors: Meeta Mistry, 

Approximate time: 90 minutes

## Learning Objectives
* Learning how to use MACS2 to compare peaks between two sample groups


## Differential binding 

Increasing number of ChIP-seq experiments are investigating transcription factor binding under multiple experimental conditions, for example, various treatment conditions, several distinct time points and different treatment dosage levels. Hence, identifying differential binding sites across multiple conditions is of practical importance in biological and medical research. 

There are various methods/tools available when investigating narrow peaks, and the choice of tool will depend heavily on your experimental design. 

![diffbind](../img/diff-peaks.png)

In our case, we are interested in identifying differences in binding between two transcription factors. For each group we have two replicates, and it would be best to use tools that make use of these replicates (i.e [DiffBind](http://bioconductor.org/packages/release/bioc/html/DiffBind.html), [ChIPComp](https://www.bioconductor.org/packages/3.3/bioc/html/ChIPComp.html)) to compute statistics reflecting how significant the changes are. 

However, in the interest of time we are going to explore tools that we have already become familiar with in this lesson. **MACS2 has a sub-command named `bdgdiff`** which can be useful for differential peak detection in the case when you don't have replicates.


## `bdgdiff` to compare peaks

This sub-command takes as input the bedgraph files that were generated as part of our peak calling. Differential peak detection is performed based on paired four bedgraph files: 

* `--t1`: MACS pileup bedGraph for condition 1 (Nanog)
* `--t2`: MACS pileup bedGraph for condition 2 (Pou5f1)
* `--c1`: MACS control lambda bedGraph for condition 1 (Nanog)
* `--c2`: MACS control lambda bedGraph for condition 2 (Pou5f1)

There are a few other parameters pertaining to the output files that we also need to specify:

* `--outdir`: the output directory
* `-o`: names for the output files 

Other parameters that can be changed (but we will leave as defaults) include:

* `-C`: logLR cutoff (default 3; likelihood ratio=1000)
* `-l`: Minimum length of differential region. Try bigger value to remove small regions (default 200)
* `-g`: Maximum gap to merge nearby differential regions. Consider a wider gap for broad marks. (default 100)
* `--d1, --d2`: Sequencing depth (# of non-redundant reads in million) for each sample

### Setting up 

Login to Orchestra and start up an interactive session:

	$ bsub -Is -q interactive bash

Navigate to the `chipseq` directory we have been working in:

	$ cd ~/ngs_course/chipseq
	
Make a directory for the ouput generated from `bdgdiff`:

	$ mkdir macs2bdgdiff
	
And finally, we will need to load the MACS2 module in order to use the sub-command:

	$ module load seq/macs/2.1.0
	
We are now ready to run the `bdgdiff` command:

```
macs2 bdgdiff --t1 macs2/Nanog-rep1_treat_pileup.bdg \
--t2 macs2/Pou5f1-rep1_treat_pileup.bdg \
--c1 macs2/Nanog-rep1_control_lambda.bdg \
--c2 macs2/Pou5f1-rep1_control_lambda.bdg \
--outdir macs2bdgdiff \
-o unique_Nanog_rep1.bed unique_Pou5f1_rep1.bed Nanog_Pou5f1_rep1common.bed 
``` 


