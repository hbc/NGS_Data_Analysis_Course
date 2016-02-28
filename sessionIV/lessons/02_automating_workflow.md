---
title: "Automating an RNA-Seq workflow"
author: "Bob Freeman, Meeta Mistry, Radhika Khetani"
date: "Wednesday, October 7, 2015"
---

## Learning Objectives:

* Automate the whole RNA-Seq workflow using a shell script
* Learn commands that make it even more flexible

## Automating the full workflow!

The easiest way to repeat the process of getting from fastqc to getting a count matrix is to capture the steps that
we've performed in a bash script. We already have 2 separate scripts for parts of this workflow, one that takes us from fastqc through trimming and a post-trimming fastqc run, and a second one that for running STAR. In this module we are going to make a new script that combines all the steps including featureCounts.

### Granting our workflow even more flexibility with a couple of new commands

We already have a good understanding of positional parameters, vectors within shell scripts, and understand the value of commenting liberally for the sake of your future self (and others who might use your script). Today, we will be putting the workflow together using 2 new commands: 

**1.** `basename`

This command will remove the full path of the file and just leave behind the filename, thus making our script more versatile. In addition, this command can make file names shorter. Let's try this out:

	$ basename ~/ngs_course/rnaseq/data/trimmed_fastq/Mov10_oe_1.subset.fq.qualtrim25.minlen35.fq
	
What if we only wanted it to return the name of the file, but without the extension .fq?

	$ basename ~/ngs_course/rnaseq/data/trimmed_fastq/Mov10_oe_1.subset.fq.qualtrim25.minlen35.fq .fq
	
What if we only wanted it to return the name of the sample?

	$ basename ~/ngs_course/rnaseq/data/trimmed_fastq/Mov10_oe_1.subset.fq.qualtrim25.minlen35.fq .fq.qualtrim25.minlen35.fq	
	
If you wanted to store the output of this command in a variable, you can write it as follows:

	$ base=$(basename ~/ngs_course/rnaseq/data/trimmed_fastq/Mov10_oe_1.subset.fq.qualtrim25.minlen35.fq .fq)
	$ echo $base

**2.** `set`

This command is essentially a debugging tool (`set -x`) that will display the command before executing it. In case of an issue with the commands in the shell script, this type of debugging lets you quickly pinpoint the step that is throwing an error. This is useful in the case where the tool is not explicitly stated in the error message, or if the error message is unclear about which tool it was created by. 

### Granting our Workflow even More Flexibility

Several changes need to be made to the last 2 scripts we made used for trimming and alignment respectively, so let's start by writing a new script with excerpts from the older ones. 

	$ cat ~/ngs_course/rnaseq/data/trimmomatic_mov10.lsf
	$ cat ~/ngs_course/rnaseq/data/trimmed_fastq/star_analysis_on_input_file.sh
	
We want to save the new script in a new directory called scripts.
	
	$ cd ~/ngs_course/rnaseq/

	$ mkdir scripts
	
	$ cd scripts

I find it easier to write a longer script in a text editor on my computer, and I suggest you do the same for this session.

```
#!bin/bash

# USAGE: sh rnaseq_analysis_on_input_file.sh <fastq files> <trimming MINLEN> <trimming TRAILING quality threshold> <number of cores>
# This script will take the location and name of a fastq file and perform the following steps on it. 
	## starting with fastqc, 
	## followed by trimming with Trimmomatic, 
	## splice-aware alignment with STAR, 
	## generation of counts associated with genes using featureCounts.

# debugging with set -x [OPTIONAL]

# set -x
```

```
# assign the command line input to new variables
fq=$1
minlen=$2
trailing=$3
cores=$4

# shorten the name of the file

fname=$(basename $fq .fq)

echo "****Running rnaseq analysis on $fname****"
```

```
# Loading all the modules and adding featureCounts to the PATH

module load seq/Trimmomatic/0.33
module load seq/STAR/2.4.0j
module load seq/fastqc/0.11.3

export PATH=/opt/bcbio/local/bin:$PATH

# make all of our output directories
	## The -p option means mkdir will create the whole path if it does not exist, and refrain from complaining if it does exist

mkdir -p ~/ngs_course/rnaseq/new_analysis
mkdir -p ~/ngs_course/rnaseq/new_analysis/trimmed_fastq
mkdir -p ~/ngs_course/rnaseq/new_analysis/STAR_alignment
mkdir -p ~/ngs_course/rnaseq/new_analysis/counts
```

```
# define variables to make modifications easier and commands shorter

trim_out=~/ngs_course/rnaseq/new_analysis/trimmed_fastq/${fname}.qualtrim${trailing}.minlen${minlen}.fq

genome=/groups/hbctraining/ngs-data-analysis2016/rnaseq/reference_data/reference_STAR 
align_out=~/ngs_course/rnaseq/new_analysis/STAR_alignment/${fname}_

gtf=~/ngs_course/rnaseq/data/reference_data/chr1-hg19_genes.gtf
counts_input_bam=~/ngs_course/rnaseq/new_analysis/STAR_alignment/${fname}_Aligned.sortedByCoord.out.bam
counts=~/ngs_course/rnaseq/new_analysis/counts/${fname}.counts
```

```
#Trimmomatic run

echo "****Trimming $fname with minimum length $minlen and trailing bases with quality threshold of $trailing****"

java -jar /opt/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads $cores -phred33 $fq \
$trim_out ILLUMINACLIP:/opt/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:$trailing MINLEN:$minlen

# FastQC on trimmed file
echo "****FastQC on trimmed fastq****"
fastqc ~/ngs_course/rnaseq/new_analysis/trimmed_fastq/${fname}.qualtrim${trailing}.minlen${minlen}.fq
```

```
# Alignment with STAR

echo "****Running STAR alignment****"

STAR --runThreadN $cores \
--genomeDir $genome \
--readFilesIn $trim_out \
--outFileNamePrefix $align_out \
--outFilterMultimapNmax 10 \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 
```

```
# Counting reads with featureCounts

echo "****Running featureCounts****"

featureCounts -T $cores -s 2 -a $gtf -o $counts $counts_input_bam
awk '{print $1"\t"$7}' $counts > $counts.txt
```

	$ vim rnaseq_analysis_on_input_file.sh
	
	$ vim submission_loop.sh

```
#!/bin/bash

for fq in ~/ngs_course/rnaseq/data/untrimmed_fastq/*.fq
do
base=$(basename $fq .subset.fq)
bsub -q mcore -n 6 -W 1:30 -R "rusage[mem=4000]" -J rnaseq_mov10.$base -o %J.out -e %J.err "sh ~/ngs_course/rnaseq/scripts/rnaseq_analysis_on_input_file.sh $fq 35 25 6"
sleep 1
done
```

> NOTE: All job schedulers are similar, but not the same. Once you understand how one works, you can transition to another one without too much trouble. They all have their pros and cons that the system administrators for your setup have taken into consideration and picked one that fits the needs of the users best. 

What you should see on the output of your screen would be the jobIDs that are returned
from the scheduler for each of the jobs that your script submitted.

You can see their progress by using the `bjobs` command (though there is a lag of
about 60 seconds between what is happening and what is reported).

Don't forget about the `bkill` command, should something go wrong and you need to
cancel your jobs.

Once all your jobs are completed, you can merge all the counts files using `paste` and `awk` together:

	$ paste ../new_analysis/counts/*.txt | head
	$ paste ../new_analysis/counts/*.txt | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12"\t"$14"\t"$16}' | head
	
	$ paste ../new_analysis/counts/*.txt | awk '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12"\t"$14"\t"$16}' > ../new_analysis/counts/all_counts.txt

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
