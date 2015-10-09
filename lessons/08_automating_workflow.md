---
title: "Automating an RNA-Seq workflow"
author: "Bob Freeman, Meeta Mistry, Radhika Khetani"
date: "Wednesday, October 7, 2015"
---

## Learning Objectives:

* Automate a workflow by grouping a series of sequential commands into a script
* Modify and submit the workflow script to the cluster


### From Sequence reads to Count matrix

That's a lot of work, and you still have five more FASTQ files to go...

- You could try running this workflow on each FASTQ file and then try and combine it all in the end. What would you have to do change in order to get this workflow to work?
- Remembering what commands *and* what parameters to type can be pretty daunting. What can
you do to help yourself out in this regard?
- How do you make sure that you are running every single file with the exact same parametrs, software versions?
- How do you keep track of all the versions and methods? (lab notebook)
- If you were to automate this process, what additional bits of information might you need?


#### Automating this Workflow with a Bash Script

The easiest way for you to be able to repeat this process is to capture the steps that
you've performed in a bash script. And you've already learned how to do this in previous
lessons. So here's a challenge...

- Using the nano text editor, create a script file that will repeat these commands
for you. You can use your command history to retrieve the commands for each step. Name your script `rnaseq_analysis_on_file.sh`. Delete your results 
directories, and run your script. Do you get all the proper output files?

- One additional command we can put in the top of the script to allow you to see what
is going on is the `set -x` bash command. This debugging tool will print out every
step before it is executed. _Insert the debugging command in your script and re-run it. How is the output different? If you're comfortable with how this looks and runs, then comment out this line._

- In order run this workflow on another file, you'll need to make changes. Copy this file,
giving the file a similar name such as `rnaseq_analysis_on_file2.sh`, and make appropriate changes to run on a control
FASTQ file called `Irrel_kd_1_qualtrim25.minlen35.fq`. _What did you have to do differently in order to get this workflow to work?_

- Reviewing your two scripts, *are there additional commonalities across scripts
or within scripts that we could optimize?*

> NOTE: Best practice would be to run the above script as a job on the scheduler. Then it is not in an interactive session, and you can set it up and exit out of the terminal. We'll do this in a few minutes.

#### Granting our Workflow More Flexibility

A couple of changes need to be made to make this script more friendly to both changes
in the workflow and changes in files. 

***More variables***

**The first major change is allowing flexibility in the input fastq file.** Thus at the start of 
the script let's capture an input parameter that must be supplied with the script name.
This input parameter will be the name of the file we want to work on:
```
    fq=$1
```
And we'll add a shortcut to where the common files are stored, for example the locations of the genome indices and the annotation file (Note that we are using the absolute paths or a path relative to our home directories, as denoted by the "~"):
```
    # location to genome reference FASTA file
    genome=/groups/hbctraining/unix_oct2015_other/reference_STAR/
    gtf=~/unix_oct2015/rnaseq_project/data/reference_data/chr1-hg19_genes.gtf
```
Make sure you are loading all the correct modules/tools for the script to run:
```    
    # set up our software environment...
    module load seq/samtools
    module load seq/htseq
```
We'll keep the output paths creation, as it looks fine. However, we will add the `-p` option this will make sure that `mkdir` will create the whole path if it does not exist, and more importantly it won't throw an error if it does exist.
```
    # make all of our output directories
    # The -p option means mkdir will create the whole path if it 
    # does not exist and refrain from complaining if it does exist
    mkdir -p ~/unix_oct2015/rnaseq_project/results/STAR
    mkdir -p ~/unix_oct2015/rnaseq_project/results/counts
```

In the script, it is a good idea to use echo for debugging/reporting to the screen
```
    echo "Processing file $fq ..."
```
We also need to use one special trick, to extract the base name of the file
(without the path and .fastq extension). We'll assign it
to the $base variable
```
    # grab base of filename for future naming
    base=$(basename $fq .qualtrim25.minlen35.fq)
    echo "basename is $base"
```
Since we've already created our output directories, we can now specify all of our
output files in their proper locations. We will assign various file names to
 variables both for convenience but also to make it easier to see what 
is going on in the command below.
```
    # set up output filenames and locations
    align_out=~/unix_oct2015/rnaseq_project/results/STAR/${base}_
    align_in=~/unix_oct2015/rnaseq_project/results/STAR/${base}_Aligned.sortedByCoord.out.bam
    counts=~/unix_oct2015/rnaseq_project/results/counts/${base}.counts
```
Our data are now staged.  We now need to change the series of command below
to use our variables so that it will run with more flexibility the steps of the 
analytical workflow
```
    # Run STAR
    STAR --runThreadN 6 --genomeDir $genome --readFilesIn $fq --outFileNamePrefix $align_out --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI NM MD AS

    # Create BAM index
    samtools index $align_in

    # Count mapped reads
    htseq-count --stranded reverse --format bam $align_in $gtf  >  $counts
```

Once you save this new script, it is ready for running:
```
$ chmod u+rwx rnaseq_analysis_on_allfiles.sh      # make it executable

$ sh rnaseq_analysis_on_allfiles.sh <name of fastq>
```

#### Running our script iteratively as a job submission to the LSF scheduler

**The above script will run in an interactive session for one file at a time. If we wanted to run this script as a job submission to LSF, and with only one command have LSF run through the analysis for all your input fastq files?**

To run the above script iteratively for all of the files on a worker node via the job scheduler, we need to create a new submission script that will need 2 important components:

* our **LSF directives** at the **beginning** of the script. This is so that the scheduler knows what resources we need in order to run our job on the compute node(s).

* a for loop that iterates through and runs the above script for all the fastq files.

Let's create a new file with nano and call it `rnaseq_analysis_on_allfiles.lsf`:
```
$ nano rnaseq_analysis_on_allfiles.lsf
```

The top of the file should look like with the LSF directives:

```
    #!/bin/bash
    #
    #BSUB -q priority		# Partition to submit to (comma separated)
    #BSUB -n 6                  # Number of cores, since we are running the STAR command with 6 threads
    #BSUB -W 1:30               # Runtime in D-HH:MM (or use minutes)
    #BSUB -R "rusage[mem=4000]"    # Memory in MB
    #BSUB -J rnaseq_mov10         # Job name
    #BSUB -o %J.out       # File to which standard out will be written
    #BSUB -e %J.err       # File to which standard err will be written
```
Then the *"for" loop*, which will make sure that all of our trimmed fastq files are being used as input:
```
    for fq in ~/unix_oct2015/rnaseq_project/data/trimmed_fastq/*.fq
    do
      rnaseq_analysis_on_allfiles.sh $fq
    done
```
We have added a few more commands at the bottom for counting; they will be executed after all the files have been processed *serially* through the `rnaseq_analysis_on_allfiles.sh` script:
```
    countmatrix=results/counts/Mov10_rnaseq_counts.txt

    # Concatenate all count files into a single count matrix
    paste results/counts/Mov10_oe_1.counts results/counts/Mov10_oe_2.counts results/counts/Mov10_oe_3.counts results/counts/Irrel_kd_1.counts results/counts/Irrel_kd_2.counts results/counts/Irrel_kd_3.counts | awk '{print$1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12}' | grep -v "^__" > $countmatrix
```

Our submission script is now complete, so save and exit out of nano. We are now ready to submit it to the queue as follows:

```
$ bsub < rnaseq_analysis_on_allfiles.lsf
```

> Now we have a count matrix for our dataset, the only thing we are missing is a header to indicate which columns correspond to which sample. We can add that in by creating a file with the header information in it:
>
>    `$ nano header.txt`
>
> Type in the following with tab separators "ID OE.1 OE.2 OE.3 IR.1 IR.2 IR.3"
>
> Now join the header to the file:
>
>    `$ cat header.txt Mov10_rnaseq_counts.txt > Mov10_rnaseq_counts_complete.txt`
>

#### Parallelizing workflow for efficiency

**The above script will run through the analysis for all your input fastq files, but it will do so in serial. We can set it up so that the pipeline is working on all the trimmed data in parallel (at the same time). This will save us a lot of time when we have realistic datasets.**

Let's make a modified version of the above script to parallelize our analysis. To do this need to modify one major aspect which will enable us to work with some of the contrainsts that this scheduler (LSF) has. We will be using a for loop for submission and putting the directives for each submission in the bsub command.

Let's make a new file called `rnaseq_analysis_on_allfiles-for_lsf.sh`

```
$ nano rnaseq_analysis_on_allfiles_for-lsf.sh
```

This file will loop through the same files as in the previous script, but the command it submits will be the actual bsub command:

```
#! /bin/bash

    for fq in ~/unix_oct2015/raw_fastq/*.fq
    do
      bsub -q priority -n 6 -W 1:30 -R "rusage[mem=4000]" -J rnaseq_mov10 -o %J.out -e %J.err "sh ~/rnaseq_analysis_on_allfiles.sh $fq"
      sleep 1
    done
```

**In the above for loop please note that after the bsub directives the `sh ~/rnaseq_analysis_on_allfiles.sh $fq` command is in quotations!**

> NOTE: All job schedulers are similar, but not the same. Once you understand how one works, you can transition to another one without too much trouble. They all have their pros and cons that the system administrators for your setup have taken into consideration and picked one that fits the needs of the users best. 

What you should see on the output of your screen would be the jobIDs that are returned
from the scheduler for each of the jobs that your script submitted.

You can see their progress by using the `bjobs` command (though there is a lag of
about 60 seconds between what is happening and what is reported).

Don't forget about the `bkill` command, should something go wrong and you need to
cancel your jobs.

#### Exercise
* Modify the `rnaseq_analysis_on_allfiles.sh` script so that one can include an additional variable to point to a results directory.
* Modify the `rnaseq_analysis_on_allfiles_for-lsf` script to work with the modified `rnaseq_analysis_on_allfiles.sh` script.

---
*The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*

