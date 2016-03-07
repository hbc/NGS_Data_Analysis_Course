---
title: "Peak calling with MACS2"
author: "Meeta Mistry"
date: "Sunday, March 6th, 2016"
---

Contributors: Meeta Mistry, Radhika Khetani

Approximate time: 90 minutes

## Learning Objectives

* Learning how to use MACS2 for peak calling
* Understanding different components of the MACS2 algorithm
* Interpretation of results from the MACS2 peak caller

## MACS2

Another popular tool for identifying transcript factor binding sites is named [Model-based Analysis of ChIP-Seq (MACS)](https://github.com/taoliu/MACS). The [MACS algorithm](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137) captures the influence of genome complexity to evaluate the significance of enriched ChIP regions. MACS improves the spatial resolution of binding sites through combining the information of both sequencing tag position and orientation. MACS can be easily used for ChIP-Seq data alone, or with control sample with the increase of specificity.


### Modeling the shift size
As we showed previously, the tag density around a true binding site should show a bimodal enrichment pattern. MACS takes advantage of this bimodal pattern to empirically model the shifting size to better locate the precise binding sites.

![model](../img/model.png)

Given a sonication size (`bandwidth`) and a high-confidence fold-enrichment (`mfold`), MACS slides 2bandwidth windows across the genome to find regions with tags more than mfold enriched relative to a random tag genome distribution. MACS randomly samples 1,000 of these high-quality peaks, separates their Watson and Crick tags, and aligns them by the midpoint between their Watson and Crick tag centers. The distance between the modes of the Watson and Crick peaks in the alignment is defined as 'd', and MACS shifts all the tags by d/2 toward the 3' ends to the most likely protein-DNA interaction sites.


### Peak detection

For experiments in which sequence depth differs between input and treatment samples, MACS linearly scales the total control tag count to be the same as the total ChIP tag count. Also, MACS allows each genomic position to contain no more than one tag and removes all the redundancies. 

To model the background noise, MACS uses a dynamic local Poisson distribution in which the lambda parameter is deduced by taking the maximum value λlocal = max(λBG, [λ1k,] λ5k, λ10k). Fold enrichment is then computed as the density of tags in a given peak compared to background λlocal parameter. This will help isolate the signal from the noise.


## Running MACS2

We will be using the newest version of this tool, MACS2. The underlying algorithm for peak calling remains the same as before, but it comes with some enhancements in functionality. One major difference on the computational end, is that for the older version (MACS) the FDR for peaks was empirically determined whereas with MACS2 the FDR is computed using the Benjamini-Hochberg method. *Note: Relaxing the q-value does not behave as expected in this case since it is partially tied to peak widths. Ideally, if you relaxed the thresholds, you would simply get more peaks but with MACS2 relaxing thresholds also results in wider peaks.* 


### Setting up

To run MACS2, we will first need to load the module:

	$ module load seq/macs/2.1.0
	

We will also need to create a directory for the output generated from MACS2:

	$ mkdir ~/ngs_course/chipseq/results/macs2
	

### MACS2 parameters

There are seven [major functions](https://github.com/taoliu/MACS#usage-of-macs2) available in MACS2 serving as sub-commands. We will only cover `callpeak` in this lesson, but if you can use `macs2 COMMAND -h` to find out more, if you are interested.

`callpeak` is the main function in MACS2 and can be invoked by typing `macs2 callpeak`. If you type this command without parameters, you will see a full description of commandline options. Here is a shorter list of the commonly used ones: 

**Input file options**

* `-t`: This is the only REQUIRED parameter for MACS
* `-c`: The control or mock data file
* `-f`: format of input file; Default is "AUTO" which will allow MACS to decide the format automatically.
* `-g`: mappable genome size which is defined as the genome size which can be sequenced; some precompiled parameters provided.

**Output arguments**

* `--outdir`: MACS2 will save all output files into speficied folder for this option
* `-n`: The prefix string for output files
* `-B/--bdg`: store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files

**Shifting model arguments**

* `-s`: size of sequencing tags. Default, MACS will use the first 10 sequences from your input treatment file to determine it
* `--bw`: The bandwidth which is used to scan the genome ONLY for model building. Can be set to the expected sonication fragment size.
* `--mfold`: upper and lower limit for model building

**Peak calling arguments**

* `-q`: q-value (minimum FDR) cutoff
* `-p`: p-value cutoff (instead of q-value cutoff)
* `--nolambda`: do not consider the local bias/lambda at peak candidate regions
* `--broad`: broad peak calling


Now that we have a feel for the different ways we can tweak our command, let's set up the command for our run on Nanog-rep1:

```
$ macs2 callpeak -t bowtie2/H1hesc_Nanog_Rep1_chr12_aln.bam \
	-c bowtie2/H1hesc_Input_Rep1_chr12_aln.bam \
 	-f BAM -g 1.3e+8 \
	--bdg --outdir macs2 \
	-n Nanog-rep1
```

The tool is quite verbose so you should see lines of text being printed to the terminal, describing each step that is being carried out. If that runs successfully, go ahead and **run the same command on the remaining samples**:

	 $ macs2 callpeak -t bowtie2/H1hesc_Nanog_Rep2_chr12_aln.bam -c bowtie2/H1hesc_Input_Rep2_chr12_aln.bam -f BAM -g 1.3e+8 --bdg --outdir macs2 -n Nanog-rep2
	 
	 $ macs2 callpeak -t bowtie2/H1hesc_Pou5f1_Rep1_chr12_aln.bam -c bowtie2/H1hesc_Input_Rep1_chr12_aln.bam -f BAM -g 1.3e+8 --bdg --outdir macs2 -n Pou5f1-rep1
	 
	 $ macs2 callpeak -t bowtie2/H1hesc_Pou5f1_Rep2_chr12_aln.bam -c bowtie2/H1hesc_Input_Rep2_chr12_aln.bam -f BAM -g 1.3e+8 --bdg --outdir macs2 -n Pou5f1-rep2

## MACS2 Output files

Let us take a look at the output files generated by `MACS2`:

	$ cd macs2/
	
	$ ls -lh
	
There should be a total of 6 files output to the results directory for each of the 4 samples:

* `_peaks.narrowPeak`: BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue
* `_peaks.xls`: a tabular file which contains information about called peaks. Additional information includes pileup and fold enrichment
* `_summits.bed`: peak summits locations for every peak. To find the motifs at the binding sites, this file is recommended
* `_model.R`: an R script which you can use to produce a PDF image about the model based on your data and cross-correlation plot
* `_control_lambda.bdg`: bedGraph format for input sample
* `_treat_pileup.bdg`: bedGraph format for treatment sample

Let's first obtain a summary of how many peaks were called in each sample. We can do this by counting the lines in the `.narrowPeak` files:

	$ wc -l macs2/*.narrowPeak

*How do these numbers compare to those generated by SPP?*

We can also take a look at the plots, but first we will have to generate it. We can use the `Rscript` command to do this:

	$ Rscript Nanog-rep1_model.r
	
Now you should see a pdf file in your current directory by the same name. Create the plots for each of the samples and move them over to your laptop using `Filezilla`. 

Open up the pdf file for Nanog-rep1. The first plot illustrates the distance between the modes from which the shift size was determined. The second plot is the cross-correlation plot, which looks different from the one generated using SPP (lower magnitude of correlation, larger shift size). This is likely due to the differences in bin size and also ranges of size shift that were evaulated.

<img src="../img/model-macs.png" width="400">




***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
