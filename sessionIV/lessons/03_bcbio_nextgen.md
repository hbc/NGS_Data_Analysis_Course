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

It is availble for installation on most Linux systems (compute clusters), and also has instructions for setup on the Cloud. It is currently installed on on the Orchestra cluster, and so we will demonstrate `bcbio-nextgen` for RNA-seq data using our Mov10 dataset as input.

The figure below describes the input (yellow), workflow (green) and output (purple) components of `bcbio`:

![bcbio](../img/bcbio-pipeline.png) 

As we work through this lesson we will introduce each component in more detail.


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

![bcbio-input](../img/bcbio-input.png) 

The files we will use as input are the **raw untrimmed FASTQ files**. Rather than than moving them and having two copies unneccesarily we can create a *symbolic link* to them. Symbolic links refer to a symbolic path indicating the abstract location of another file. The syntax is as follows:

	ln -s {/path/to/file-name} {link-name}

In our situation we would use:

	ln -s ~/ngs_course/rnaseq/data/untrimmed_fastq/ raw-data

This would allow us to access the data directly from another directory without having to provide paths.

In addition to the data files, `bcbio` requires a **comma separated value file containing sample metadata**. The first column must contain the header `samplename` which corresponds to the FASTQ filenames you are running the analysis on. You can add a `description` column to change the sample name originally supplied by the file name, to this value (i.e. a short name). And finally, any columns that follow can contain additional information on each sample.

We have created this file for you, you will need to copy it over to your current directory.

	cp /groups/hbctraining/ngs-data-analysis2016/rnaseq/bcbio-rnaseq/mov10_project.csv .
	
Each line in the file corresponds to a sample, and each column has information about the samples.

```
samplename,description,condition
Irrel_kd_1.subset.fq,Irrel_kd_1,control
Irrel_kd_2.subset.fq,Irrel_kd_2,control
Irrel_kd_3.subset.fq,Irrel_kd_3,control
Mov10_oe_1.subset.fq,Mov10_oe_1,overexpression
Mov10_oe_2.subset.fq,Mov10_oe_2,overexpression
Mov10_oe_3.subset.fq,Mov10_oe_3,overexpression
```

The final requirement is a **configuration template**, which will contain details on the analysis options. The template file is used in combination with the metadata file, and the FASTQ files to create a **config file** which is ultimately the input for `bcbio`.

You can start with one of the provided [best-practice templates](https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates) and modify it as required, or you can create your own. We have created a template for you based on the experimental details. Copy it over and then use `less` to take a look at what is inside.

	cp /groups/hbctraining/ngs-data-analysis2016/rnaseq/bcbio-rnaseq/mov10-template.yaml
	less mov10-template.yaml
	
```
# Template for human RNA-seq using Illumina prepared samples
---
details:
  - analysis: RNA-seq
    genome_build: hg19
    algorithm:
      aligner: star
      quality_format: standard
      trim_reads: read_through
      adapters: [truseq, polya]
      strandedness: firststrand 
upload:
  dir: ../final
```

You should observe indentation which is characteristic of the **YAML file format** (YAML Ain't Markup Language). YAML is a human friendly data serialization standard for all programming languages. It takes concepts from languages such as C, Perl, and Python and ideas from XML. The data structure hierarchy is maintained by **outline indentation**.

The configuration template defines `details` of each sample to process, including : `analysis` (i.e. RNA-seq, chipseq, variant), `genome_build`, and `algorithm` specifics for each tool that is being used in the workflow. Other details include `metadata`, which will be added via the `.csv` when creating the final config file. At the end of the template we define `upload` which is for the final ouput from `bcbio`. To find out more on the details that can be added to your YAML, check out the [readthedocs](https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#sample-information). 	
	 
We can now apply this template to all samples in our dataset. To do this we use the	template workflow command, which takes in the template, the metadata and the samples:  

	bcbio_nextgen.py -w template mov10-template.yaml mov10_project.csv raw_fastq/*.fq
	
Upon completion of the command you should see the following output:

```
Template configuration file created at: /home/mm573/ngs_course/rnaseq/bcbio-rnaseq/mov10_project/config/mov10_project-template.yaml
```

If you take a look in your current directory, you will also find that a new directory has been created by the same name as your csv file `mov10_project`. Inside that directory you will find the following directory structure:

```
mov10_project/
├── config
└── work
```

## `bcbio`: workflow




***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
