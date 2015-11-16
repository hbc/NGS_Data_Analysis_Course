---
title: "Getting your project started"
author: "Jason Williams, Bob Freeman, Meeta Mistry"
date: "Wednesday, October 7, 2015"
---

Approximate time: 30 minutes

## Learning Objectives

* Understand the experiment and its objectives
* Setting up for a NGS data analysis 
* Learning best practices for project organization


## Understanding the dataset
The dataset we are using is part of a larger study described in [Kenny PJ et al, Cell Rep 2014](http://www.ncbi.nlm.nih.gov/pubmed/25464849). The authors are investigating interactions between various genes involved in Fragile X syndrome, a disease in which there is aberrant production of the FMRP protein. FMRP has been linked to the microRNA pathway, as it has been shown to be involved in miRNA mediated translational suppresion. **The authors sought to show that FMRP associates with the RNA helicase MOV10, that is also associated with the microRNA pathway.**

From this study we are using the [RNA-Seq](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50499) data which is publicly available in the [SRA](http://www.ncbi.nlm.nih.gov/sra). The RNA was extracted from HEK293F cells that were transfected with a MOV10 transgene and normal control cells. Using this data, we will evaluate transcriptional patterns associated with MOV10 overexpression. The libraries for this dataset are stranded and were generated using the dUTP method. Sequencing was carried out on the Illumina HiSeq-2500 for 100bp single end reads. The full dataset was sequenced to ~40 million reads per sample, but for this workshop we will be looking at a small subset on chr1 (~300,000 reads/sample). For each group we have three replicates as described in the figure below.


![Automation](../img/exp_design.png)

 
## Getting your project started

Project organization is one of the most important parts of a sequencing project, but is often overlooked in the excitement to get a first look at new data. While it's best to get yourself organized before you begin analysis, it's never too late to start.

You should approach your sequencing project in a very similar way to how you do a biological experiment, and ideally, begins with experimental design. **We're going to assume that you've already designed a beautiful sequencing experiment to address your biological question, collected appropriate samples, and that you have enough statistical power.** We will cover experimental design in more detail as they apply to the different workflows. 

For all of those steps taken in the wetlab (collecting specimens, extracting DNA, prepping your samples) you've likely kept a lab notebook that details how and why you did each step, but documentation doesn't stop at the sequencer! 

Every computational analysis you do is going to spawn many files, and inevitability, you'll want to run some of those analysis again. **Genomics projects can quickly accumulate hundreds of files across tens of folders.** Do you remember what PCR conditions you used to create your sequencing library? Probably not. Similarly, you probably won't remember whether your best alignment results were in Analysis1, AnalysisRedone, or AnalysisRedone2; or which quality cutoff you used.

Luckily, recording your computational experiments is even easier than recording lab data. Sensible file names will make your analysis traversable by you and your collaborators, and writing the methods section for your next paper will be a breeze. Let's look at the best practices for documenting your genomics project. 

Your future self will thank you.

> Adapted from [Data Organization Best Practices](https://github.com/datacarpentry/organization-genomics/blob/gh-pages/GoodBetterBest.md)

### Setting up the filesystem

In this next exercise we will setup a filesystem for the project we will be using over the next few days. \

Inspired by the guide below, we will start by creating a directory that we can use for the rest of the workshop:

First, make sure that you are in your home directory,

```
$ pwd
```
this should give the result: `/home/user_name`

* **Tip** Remember, when we give a command, rather than copying and pasting, just type it out. Also the '$' indicates we are at the command prompt, do not include that in your command. 
* **Tip** If you were not in your home directory, the easiest way to get there is to enter the command *cd* - which always returns you to home. 

Now, make a directory for your project within the `unix_oct2015` folder using the `mkdir` command

```
$ mkdir unix_oct2015/rnaseq_project
```

Next you want to set up the following directory structure within your project directory to keep files organized. 

```
rnaseq_project/
├── data
├── meta
├── results
└── docs

```
You can do this by changing into `rnaseq_project` and then using `mkdir` to create the four directories.

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
if you have created these directories, you should get the following output from that command:

```
/data  /docs  /meta  /results

```



### Best practices for naming files


### Document your activity on the project

> It is also useful to have README file within your project directory. This file will usually contain a quick one line summary about the project and any other lines that follow will describe the files/directories found within it. 
> 
> Take a moment to create a README for `rnaseq_project` (hint: use nano to create the file). Give a short description of the project and brief descriptions of the types of file you would be storing within each of the sub-directories.

The `history` command is a convenient way to document the all the commands you have used while analyzing and manipulating your project. Let's document the work we have done to create these folders. 

View the commands that you have used so far during this session using history:

```
$ history
```

The history likely contains many more commands that you have used just for these projects. Let's view the last several commands so that focus on just what we need for the project. 

 View the last n lines of your history (where n = approximately the last few lines you think relevant - for our example we will use the last 7:

```
$ history | tail -n7
```

As you may remember from the shell lesson, the pipe '|' sends the output of history to the next program, in this case, tail. We have used the -n option to give the last 7 lines.

Using your knowledge of the shell use the append redirect `'>>'` to create a file called **unix_workshop_log_XXXX_XX_XX.txt** (Use the four-digit year, two-digit month, and two digit day, e.g. unix_workshop_log_2015_10_08.txt)


You may have noticed that your history may contain the *history* command itself. To remove this redundancy from our log, lets use the *nano* text editor to fix the file:

```
$ nano unix_workshop_log_
```

From the nano screen, you should be able to use your cursor to navigate, type, and delete any redundant lines. 

Add a dateline and comment to the line where you have created the directory e.g. 

```
# 2015_10_08 
```

```
# Created sample directories for the Intro to Unix workshop
```

6. Next, remove any lines of the history that are not relevant. Just navigate to those lines and use your delete key. 
7. Close nano by hitting 'Control' and the 'X' key at the same time; notice in nano this is abbreviated '\^X'; nano will ask if you want to save; hit 'Y' for yes. When prompted for the 'File Name to Write' we can hit 'Enter' to keep the same name and save. 
8. Now that you have created the file, move the file to 'rnaseq_project/docs'



###References
[A Quick Guide to Organizing Computational Biology Projects] (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)


---
*The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*

