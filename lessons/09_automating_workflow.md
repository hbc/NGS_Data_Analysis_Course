# Lesson

Automating an RNA-Seq workflow
===================

Learning Objectives:
-------------------
### What's the goal for this lesson?

* Automate a workflow by grouping a series of sequential commands into a script
* Modify and submit the workflow script to the cluster



### From Sequence reads to Count matrix

That's a lot of work, yes? And you still have five more FASTQ files to go...

- You could try running this workflow on each FASTQ file andthen try and combine it all in the end. What would you have to do change
in order to get this workflow to work?
- Remembering what commands *and* what parameters to type can be pretty daunting. What can
you do to help yourself out in this regard?
- If you were to automate this process, what additional bits of information might you need?


#### Automating this Workflow with a Bash Script

The easiest way for you to be able to repeat this process is to capture the steps that
you've performed in a bash script. And you've already learned how to do this in previous
lessons. So...

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


#### Granting our Workflow More Flexibility

A couple of changes need to be made to make this script more friendly to both changes
in the workflow and changes in files. 

**The first major change is allowing a change in the filename.** Thus at the start of 
the script let's capture an input parameter that must be supplied with the script name.
This input parameter will be the name of the file we want to work on:

     fq=$1

And we'll add a shortcut to where the common files are stored, for example: the location of the trimmed FASTQ and locations of the genome indices and the annotation file:

     # location to FASTQ file
     fileIn=data/untrimmed_fastq

     # location to genome reference FASTA file
     genome=/groups/hbctraining/unix_oct2015_other/reference_STAR/
     gtf=data/reference_data/

Make sure you are loading all the correct modules/tools for the script to run:
    
    # set up our software environment...
    module load seq/samtools
    module load seq/htseq

We'll keep the output paths creation, as it looks fine. (Though really, we could
put results/ in a variable and declare that at the top, so we can change where the
results will be as well. We'll leave that for an optional exercise)

     # make all of our output directories
     mkdir -p results/STAR
     mkdir -p results/counts


In the script, it is a good idea to use echo for debugging/reporting to the screen

    echo "Processing file $fq ..."

We also need to use one special trick, to extract the base name of the file
(without the path and .fastq extension). We'll assign it
to the $base variable

    # grab base of filename for future naming
    base=$(basename $fq .qualtrim25.minlen35.fq)
    echo "basename is $base"

Since we've already created our output directories, we can now specify all of our
output files in their proper locations. We will assign various file names to
 variables both for convenience but also to make it easier to see what 
is going on in the command below.

    # set up output filenames and locations
    align_out=results/STAR/$base
    align_in=results/STAR/${base}__Aligned.sortedByCoord.out.bam
    counts=results/counts/${base}.counts

Our data are now staged.  We now need to change the series of command below
to use our variables so that it will run with more flexibility the steps of the 
analytical workflow

```
# Run STAR
STAR --runThreadN 6 --genomeDir $genome --readFilesIn ${fileIn}/${fq} --outFileNamePrefix $align_out --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI NM MD AS

# Create BAM index
samtools index $align_in

# Count mapped reads
htseq-count --stranded reverse --format bam $align_in $gtf  >  $counts

```

This new script is now ready for running:
	
	sh run_variant_call_on_file.sh <name of fastq>

#### Parallelizing workflow for efficiency

To run the same script on a worker node on the cluster via the job scheduler, we need to add our **LSF directives** at the **beginning** of the script. This is so that the scheduler knows what resources we need in order to run our job on the
compute node(s). 

Copy the `rnaseq_analysis_on_file.sh` file and give it a new name `rnaseq_analysis_on_file.lsf`. Add the LSF directives to the beginning of the file.

So the top of the file should look like:

    #!/bin/bash
    #
    #BSUB -q priority		# Partition to submit to (comma separated)
    #BSUB -n 1                  # Number of cores
    #BSUB -W 1:30               # Runtime in D-HH:MM (or use minutes)
    #BSUB -R "rusage[mem=10000]"    # Memory in MB
    #BSUB -J rnaseq_mov10         # Job name
    #BSUB -o %J.out       # File to which standard out will be written
    #BSUB -e %J.err       # File to which standard err will be written


What we'd like to do is run this script on a compute node for every trimmed FASTQ -- pleasantly parallelizing our workflow. And now is where we'll use the for loop with the power of the cluster: 

    for fq in data/trimmed_fastq/*.fastq
    do
      bsub rnaseq_analysis_on_file.sh $fq
      sleep 1
    done

What you should see on the output of your screen would be the jobIDs that are returned
from the scheduler for each of the jobs that you submitted.

You can see their progress by using the `bjobs` command (though there is a lag of
about 60 seconds between what is happening and what is reported).

Don't forget about the `bkill` command, should something go wrong and you need to
cancel your jobs.

#### Exercise
* Change the script so that one can include an additional variable to point to
 a results directory.
