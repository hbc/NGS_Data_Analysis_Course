---
title: "Getting your project started"
author: "Jason Williams, Bob Freeman, Meeta Mistry"
date: "Wednesday, October 7, 2015"
---

Approximate time: 30 minutes

## Learning Objectives

* Understand the experiment and its objectives
* Setting up your project space for an NGS workflow
* Learning best practices for NGS analysis


## Understanding the dataset
The dataset we are using is part of a larger study described in [Kenny PJ et al, Cell Rep 2014](http://www.ncbi.nlm.nih.gov/pubmed/25464849). The authors are investigating interactions between various genes involved in Fragile X syndrome, a disease in which there is aberrant production of the FMRP protein. FMRP has been linked to the microRNA pathway, as it has been shown to be involved in miRNA mediated translational suppresion. **The authors sought to show that FMRP associates with the RNA helicase MOV10, that is also associated with the microRNA pathway.**

From this study we are using the [RNA-Seq](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50499) data which is publicly available in the [SRA](http://www.ncbi.nlm.nih.gov/sra). The RNA was extracted from HEK293F cells that were transfected with a MOV10 transgene and normal control cells. Using this data, we will evaluate transcriptional patterns associated with MOV10 overexpression. The libraries for this dataset are stranded and were generated using the dUTP method. Sequencing was carried out on the Illumina HiSeq-2500 for 100bp single end reads. The full dataset was sequenced to ~40 million reads per sample, but for this workshop we will be looking at a small subset on chr1 (~300,000 reads/sample). For each group we have three replicates as described in the figure below.


![Automation](../img/exp_design.png)

 
## Getting your project started

Project organization is one of the most important parts of a sequencing project, but is often overlooked in the excitement to get a first look at new data. While it's best to get yourself organized before you begin analysis, it's never too late to start.

You should approach your sequencing project in a very similar way to how you do a biological experiment, and ideally, begins with experimental design. **We're going to assume that you've already designed a beautiful sequencing experiment to address your biological question, collected appropriate samples, and that you have enough statistical power.** We will cover experimental design in more detail as they apply to the different workflows. 

For all of those steps taken in the wetlab (collecting specimens, extracting DNA, prepping your samples) you've likely kept a lab notebook that details how and why you did each step, but documentation doesn't stop at the sequencer! 

Every computational analysis you do is going to spawn many files, and inevitability, you'll want to run some of those analysis again. **Genomics projects can quickly accumulate hundreds of files across tens of folders.** Do you remember what PCR conditions you used to create your sequencing library? Probably not. Similarly, you probably won't remember whether your best alignment results were in Analysis1, AnalysisRedone, or AnalysisRedone2; or which quality cutoff you used.

Luckily, recording your computational experiments is even easier than recording lab data. Sensible file names will make your analysis traversable by you and your collaborators, and writing the methods section for your next paper will be a breeze. Let's look at the best practices for organizing your genomics project. 

Your future self will thank you.


### Setting up the filesystem

To prepare for the overwhelming number of files that will be generated during your NGS workflow, it is best to set up a directory structure such that you have a designated place for files when you encounter them. In this next exercise we will setup a directory structure for the project we will be using over the next few days.

First, make sure that you are in your home directory,

```
$ pwd
```
this should give the result: `/home/user_name`

**Tip** If you were not in your home directory, the easiest way to get there is to enter the command *cd* - which always returns you to home. 

Now, make a directory for your project within the `unix_oct2015` folder using the `mkdir` command

```
$ mkdir unix_oct2015/rnaseq_project
```

Next you want to set up the following structure within your project directory to keep files organized:

```
rnaseq_project/
├── data
├── meta
├── results
└── logs

```

This is a generic structure and can be tweaked based on personal preferences. A brief description of what might be contained within the different sub-directories is provided below:

* **`data/`**: This folder is usually reserved for any raw data files that you start with. For example, this is where the original FASTQ files (data you get from the sequencer) would reside. It is best practice to always have a copy of the data in a folder in it's raw form; as you will notice that usually a workflow is run a few times before we get it completely right.
* **`meta/`**: This folder contains any information that describes the samples you are using, which we often refer to as metadata. Usually this comes in the form of a tab-delimited/Excel file in which each row corresponds to a sample (listed using the filename for that sample in the raw data collection), and columns that follow would contain any other pertinent information for the sample (i.e sample class, demographic factors, sequencer specific information). An example of a metadata file is shown below:

![metadata](../img/metadata_example.png)

* **`results/`**: This folder will contain the output from the different tools you implement in your workflow. In some cases, you will simply have the results file as ouput but with other tools you will find a large number of intermediate files are generated. To stay organized, you should create sub-folders specific to each tool/step of the workflow. 
* **`logs/`**: It is important to keep track of the commands you run and the specific pararmeters you used, but also to have a record of any standard output that is generated while running the command. This will allow you to go back to your recorded logfiles to explore additional information (e.g., how many adapters were removed, how many reads did not align). Different tools have different ways of reporting log messages and you might have to experiment a bit to figure out what output to capture: you can redirect standard output with the `>` symbol which is equivalent to `1> (standard out)`; other tools might require you to use `2>` to re-direct the standard error instead. 
 

Let's create the directory structure for our  changing into `rnaseq_project` and then using `mkdir` to create the four directories.

```
$ cd unix_oct2015/rnaseq_project
$ mkdir data
$ mkdir meta
$ mkdir results
$ mkdir docs

``` 

Verify that you have created the directories:

```
$ ls -F
```
 
If you have created these directories, you should get the following output from that command:

```
/data  /docs  /meta  /results

```

### Document your activity on the project

Keeping notes on what happened in what order, what was done and by whom, is essential for reproducible research.  It is essential for good science.  If you don’t keep good notes, then you will forget what you did pretty quickly, and if you don’t know what you did, no-one else has a chance. After setting up the filesystem it is useful to have a README file within your project directory. This file will usually contain a quick one line summary about the project and any other lines that follow will describe the files/directories found within it. Within each sub-directory you can also include README files to describe the files that were generated. 
 
> Take a moment to create a `README.txt` for `rnaseq_project` (hint: use nano to create the file). Give a short description of the project with today's date and a brief descriptions of the types of file you intend to store within each of the sub-directories.


To keep track of the commands you have used while analyzing your data, the `history` command is very convenient. We haven't gotten to any data analysis just yet, but as an example we can document the commands we have used to create these folders. 

To view the commands that you have used so far during this session using history:

```
$ history
```

The history likely contains many more commands that you have used just for these projects. Let's view the last several commands so that focus on just what we need for the project. View the last n lines of your history (where n = approximately the last few lines you think relevant - for our example we will use the last 7:

```
$ history | tail -n7
```

As you may remember from the shell lesson, the pipe '|' sends the output of history to the next program, in this case, tail. We have used the -n option to give the last 7 lines. Using your knowledge of the shell use the append redirect `'>>'` to create a file called **ngs_workshop_log_XXXX_XX_XX.txt** (Use the four-digit year, two-digit month, and two digit day, e.g. ngs_workshop_log_2015_10_08.txt)


You may have noticed that your history may contain the *history* command itself. To remove this redundancy from our log, lets use the `nano` text editor to fix the file. From the nano screen, you should be able to use your cursor to navigate, type, and delete any redundant lines. 

Add a dateline and comment above the lines of history:

```
# 2015_10_08 
# Created sample directories for the Intro to Unix workshop
```

Next, remove any lines of the history that are not relevant. Just navigate to those lines and use your delete key. Close nano by hitting 'Control' and the 'X' key at the same time. Now that you have created the file, move the file to `rnaseq_project/logs`

### Naming files

A few months from now, you may not remember what you were up to when you created a particular set of files. Below is a short list of things we suggest when it comes to file naming:

1. **Keep sample names short and meaningful.** If required, include some form a long explanation for the sample names (i.e comment lines at the top of the metadata file, or add it in your README file).
2. Have **unique sample names** and try to avoid names that look like dates (Dec14), times (AM1245) and other things that Excel might auto-convert. 
3. **Remove spaces and punctuation.** When working on the command line, spaces in file names make everything exponentially more difficult. Replace all your spaces with under_scores and avoid the use of any special characters.

## Best practices for NGS Analysis 

Ok so now you are all set up to start your analyses! You have set up your space in a way such that someone unfamiliar with your project should be able to look at your computer files and understand in detail what you did and why. Now before we move on to any actual data, we have a few words of wisdom to impart upon you:


1. **Make sure to use the appropriate software.** Do your research and find out what is best for the data you are working with. Don't just work with tools that you are able to easily install. Also, make sure you are using the most up-to-date versions! If you run out-of-date software, you are probably introducing errors into your workflow; and you may be missing out on more accurate methods.

2. **Keep up with the literature.** Bioinformatics is a fast-moving field and it's always good to stay in the know about recent developments. This will help you determine what is appropriate and what is not.  

3. **Do not re-invent the wheel.** If you run into problems, more often than not someone has already encountered that same problem. A solution is either already available or someone is working on it -- so find it!

4. **Testing is essential.** If you are using a tool for the first time, test it out on a single sample or a subset of the data before running your entire dataset through. This will allow you to debug quicker and give you a chance to also get a feel for the tool and the different parameters.



###References
* [A Quick Guide to Organizing Computational Biology Projects] (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)
* [Data Organization Best Practices](https://github.com/datacarpentry/organization-genomics/blob/gh-pages/GoodBetterBest.md)
* [The five habits of bad bioinformaticians](https://biomickwatson.wordpress.com/2015/11/16/the-five-habits-of-bad-bioinformaticians/)



---
*The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*

