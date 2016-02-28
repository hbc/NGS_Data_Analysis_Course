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

# assign the command line input to new variables
fq=$1
minlen=$2
trailing=$3
cores=$4

# shorten the name of the file

fname=$(basename $fq .fq)

echo "****Running rnaseq analysis on $fname****"

# Loading all the modules and adding featureCounts to the PATH

module load seq/Trimmomatic/0.33
module load seq/STAR/2.4.0j
module load seq/fastqc/0.11.3

export PATH=/opt/bcbio/local/bin:$PATH

# make all of our output directories
# The -p option means mkdir will create the whole path if it does not exist, and refrain from complaining if it does exist

mkdir -p ~/ngs_course/rnaseq/new_analysis/

mkdir -p ~/ngs_course/rnaseq/new_analysis/trimmed_fastq

mkdir -p ~/ngs_course/rnaseq/new_analysis/STAR_alignment

mkdir -p ~/ngs_course/rnaseq/new_analysis/counts


# setup variables

trim_out=~/ngs_course/rnaseq/new_analysis/trimmed_fastq/${fname}.qualtrim${trailing}.minlen${minlen}.fq

genome=/groups/hbctraining/ngs-data-analysis2016/rnaseq/reference_data/reference_STAR 
align_out=~/ngs_course/rnaseq/new_analysis/STAR_alignment/${fname}_

gtf=~/ngs_course/rnaseq/data/reference_data/chr1-hg19_genes.gtf
counts_input_bam=~/ngs_course/rnaseq/new_analysis/STAR_alignment/${fname}_Aligned.sortedByCoord.out.bam
counts=~/ngs_course/rnaseq/new_analysis/counts/${fname}.counts

#Trimmomatic run

echo "****Trimming $fname with minimum length $minlen and trailing bases with quality threshold of $trailing****"

java -jar /opt/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads $cores -phred33 $fq \
$trim_out ILLUMINACLIP:/opt/Trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:$trailing MINLEN:$minlen

# FastQC on trimmed file
echo "****FastQC on trimmed fastq****"
fastqc ~/ngs_course/rnaseq/new_analysis/trimmed_fastq/${fname}.qualtrim${trailing}.minlen${minlen}.fq

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

# Counting reads with featureCounts

echo "****Running featureCounts****"

featureCounts -T $cores -s 2 -a $gtf -o $counts $counts_input_bam
awk '{print $1"\t"$7}' $counts > $counts.txt
```

Next, we'll initialize variables that contain the paths to where the common files are stored and then use the variable names (with a `$`) in the actual commands later in the script. This is a shortcut for when you want to use this script for a dataset that used a different genome, e.g. mouse; you'll just have to change the contents of these variable at the beginning of the script.

Let's add 2 variables named "genome" and "gtf", these will contain the locations of the genome indices and the annotation file respectively:

    # location of genome reference FASTA and index files + the gene annotation file
    genome=/groups/hbctraining/unix_oct2015_other/reference_STAR/
    gtf=~/unix_oct2015/rnaseq_project/data/reference_data/chr1-hg19_genes.gtf

Next, make sure you load all the modules for the script to run. This is important so your script can run independent of any "prep" steps that need to be run beforehand:
    
    # set up our software environment...
    module load seq/samtools
    module load seq/htseq

We'll keep the output directory creation, however, we will add the `-p` option this will make sure that `mkdir` will create the directory only if it does not exist, and it won't throw an error if it does exist.
```
    # make all of our output directories
    # The -p option means mkdir will create the whole path if it 
    # does not exist and refrain from complaining if it does exist
    mkdir -p ~/unix_oct2015/rnaseq_project/results/STAR
    mkdir -p ~/unix_oct2015/rnaseq_project/results/counts
```

In the script, it is a good idea to use echo for debugging/reporting to the screen (you can also use `set -x`):
```
    echo "Processing file $fq ..."
```
We also need to use one special trick, to extract the base name of the file
(without the path and .fastq extension). We'll assign it to the $base variable using the basename command:
```
    # grab base of filename for future naming
    base=$(basename $fq .qualtrim25.minlen35.fq)
    echo "basename is $base"
```
There are 2 new things of note above:

1. the `basename` command: this command takes a path or a name and trims away all the information before the last `\` and if you specify the string to clear away at the end, it will do that as well. In this case, if the variable `$fq` contains the path *"unix_oct2015_other/trimmed_fastq/Mov10_oe_1.qualtrim25.minlen35.fq"*, `basename $fq .qualtrim25.minlen35.fq` will output "Mov10_oe_1".
2. to assign this value to the `base` variable, we place the `basename...` command in parentheses and put a `$` outside. This syntax is necessary for assigning the output of a command to a variable.

Since we've already created our output directories, we can now specify all of our
output files in their proper locations. We will assign various file names to
 variables both for convenience but also to make it easier to see what 
is going on in the command below.
```
    # set up output filenames and locations
    align_out=~/unix_oct2015/rnaseq_project/results/STAR/${base}_
    counts_input_bam=~/unix_oct2015/rnaseq_project/results/STAR/${base}_Aligned.sortedByCoord.out.bam
    counts=~/unix_oct2015/rnaseq_project/results/counts/${base}.counts
```
Our variables are now staged. We now need to modify the series of commands starting with STAR throught to counts (htseq-count)
to use them so that it will run the steps of the analytical workflow with more flexibility:

    # Run STAR
    STAR --runThreadN 6 --genomeDir $genome --readFilesIn $fq --outFileNamePrefix $align_out --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI NM MD AS

    # Create BAM index
    samtools index $align_in

    # Count mapped reads
    htseq-count --stranded reverse --format bam $counts_input_bam $gtf  >  $counts


Once you save this new script, it is ready for running:
```
$ chmod u+rwx rnaseq_analysis_on_input_file.sh      # make it executable, this is good to do, even if your script runs fine without it to ensure that it always does and you are able to tell that it's an executable shell script.

$ sh rnaseq_analysis_on_input_file.sh <name of fastq>
```

### Parallelizing workflow for efficiency

**The above script will run through the analysis for all your input fastq files, but it will do so in serial. We can set it up so that the pipeline is working on all the trimmed data in parallel (at the same time). This will save us a lot of time when we have realistic datasets.**

Let's make a modified version of the above script to parallelize our analysis. To do this need to modify one major aspect which will enable us to work with some of the contrainsts that this scheduler (LSF) has. We will be using a for loop for submission and putting the directives for each submission in the bsub command.

Let's make a new file called `rnaseq_analysis_on_allfiles-for_lsf.sh`. Note this is a normal shell script.


	$ nano rnaseq_analysis_on_allfiles_for-lsf.sh


This file will loop through the same files as in the previous script, but the command it submits will be the actual bsub command:


	#! /bin/bash

    for fq in ~/unix_oct2015/raw_fastq/*.fq
    do
      bsub -q priority -n 6 -W 1:30 -R "rusage[mem=4000]" -J rnaseq_mov10 -o %J.out -e %J.err "sh rnaseq_analysis_on_input_file.sh $fq"
      sleep 1
    done


**In the above for loop please note that after the bsub directives the `sh rnaseq_analysis_on_input_file.sh $fq` command is in quotes!**

> NOTE: All job schedulers are similar, but not the same. Once you understand how one works, you can transition to another one without too much trouble. They all have their pros and cons that the system administrators for your setup have taken into consideration and picked one that fits the needs of the users best. 

What you should see on the output of your screen would be the jobIDs that are returned
from the scheduler for each of the jobs that your script submitted.

You can see their progress by using the `bjobs` command (though there is a lag of
about 60 seconds between what is happening and what is reported).

Don't forget about the `bkill` command, should something go wrong and you need to
cancel your jobs.

---
*The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*

