---
title: "NGS pipelines: bcbio-nextgen"
author: "Meeta Mistry"
date: "Thursday, February 25, 2015"
---

Approximate time: 


## NGS pipelines

As you can see from our RNA-seq lessons so far, the analysis workflow is a multi-step process. We learned what is involved in running each individual step, and the details on inputs and outputs. Finally, we demonstrated how to combine the different steps into a single script for automation of the entire workflow from start to finish.

An alternative to creating your own pipeline for the analysis of your next-generation sequencing data, it to use an existing one. There are a number of pipelines available both commerical and academic, with some that are specific to a particular NGS experiment (i.e [RNA-seq only](http://www.biomedcentral.com/content/pdf/1471-2164-16-S6-S3.pdf), [viral NGS analysis](http://viral-ngs.readthedocs.org/en/latest/)).

The pipeline we will be presenting here is [bcbio-nextgen](https://bcbio-nextgen.readthedocs.org/en/latest/).


## `bcbio-nextgen`

bcbio-nextgen is a shared community resource for handling the data processing of various NGS experiments including variant calling, mRNA-seq, small RNA-seq and ChIP-seq. It is an open-source toolkit created by Brad Chapman with a lot of support (testing, bug reports etc.) and development from the large [user community](https://bcbio-nextgen.readthedocs.org/en/latest/contents/introduction.html#users).

> *"A piece of of software is being sustained if people are using it, fixing it, and improving it rather than replacing it"*
> -[Software Carpentry](http://software-carpentry.org/blog/2014/08/sustainability.html)
>

bcbcio-nextgen provides *best-practice* piplelines with the goal of being:

* Scalable: Handles large datasets and sample populations on distributed heterogeneous compute environments
* Well-documented
* Easy(ish) to use: Tools come pre-configured
* Reproducible: Tracks configuration, versions, provenance and command lines
* Analyzable: Results feed into downstream tools to make it easy to query and visualize

It is availble for installation on most Linux systems (compute clusters), and instructions for setup on the Cloud. It is currently installed on on the Orchestra cluster, and so we will demonstrate `bcbio-nextgen` for RNA-seq data using our Mov10 dataset as input.

The figure below describes the input (yellow), workflow (green) and output (purple) components of `bcbio`:

![bcbio](../img/bcbio-pipeline.png) 

As we work through this lesson we will introduce each component in more detail .


## Setting up

Let's get started by logging on to Orchestra and starting an interactive session:

	bsub -Is -q interactive bash
	
The first thing we need to do in order to run `bcbio`, is setup some environment variables. Rather than just modifying them in the command-line, we will be adding it to our `.bashrc` file which is  located in your home directory. The `.bashrc` is a shell script that Bash runs whenever it is started interactively. You can put any command in that file that you could type at the command prompt, and is generally used to set an environment and customize things to your preferences.

Open up your `.bashrc` using `vim` and add in the following:

	# Environment variables for running bcbio
	export PATH=/opt/bcbio/local/bin:$PATH
	export LD_LIBRARY_PATH=/opt/bcbio/local/lib:$LD_LIBRARY_PATH
	export PERL5LIB=/opt/bcbio/local/lib/perl5:${PERL5LIB}
 
Close and save the file. Finally, let's set up the project structure. Change directories into `~/ngs_course/rnaseq` and make a directory called `bcbio-rnaseq`:

	cd ~/ngs_course/rnaseq
	mkdir bcbio-rnaseq



## `bcbio`: Inputs

There are three things required as input for your `bcbio` run:

1. FASTQ/BAM files
2. Metadata file (.csv)
3. Configuration file

![bcbio-input](../img/bcbio-input.png) 

The first of this list, we already have. The files we will use as input are the **raw untrimmed FASTQ files**. Rather than than moving them and having two copies unneccesarily we can create a *symbolic link* to them. Symbolic links refer to a symbolic path indicating the abstract location of another file. The syntax is as follows:

	ln -s {/path/to/file-name} {link-name}

In our situation we would use:

	ln -s ~/ngs_course/rnaseq/data/untrimmed_fastq/ raw-data

This would allow us to access the data directly from another directory without having to provide paths.

In addition to the data files, `bcbio` requires a **comma separated value file containing sample metadata**. The first column must contain the header `samplename` which corresponds to the FASTQ filenames you are running the analysis on. You can add a `description` column to change the sample name originally supplied by the file name, to this value (i.e. a short name). And finally, any columns that follow can contain additional information on each sample.

We have created this file for you, you will need to copy it over to your current directory.

	cp /groups/hbctraining/ngs-data-analysis2016/rnaseq/bcbio-rnaseq/mov10_project.csv .
	

The final requirement is a **configuration file**, which will contain details about the set of samples to process, including input files and analysis options. The conifguration file is provided in **YAML file format**, which stands for YAML Ain't Markup Language. YAML is a human friendly data serialization standard for all programming languages. It takes concepts from languages such as C, Perl, and Python and ideas from XML. The data structure hierarchy is maintained by outline indentation.


To create this configuration file, bcbio-nextgen provides a utility to create configuration files for multiple sample inputs using a base template. Start with one of the best-practice templates, or define your own (mov10_template.yaml).


	





***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
