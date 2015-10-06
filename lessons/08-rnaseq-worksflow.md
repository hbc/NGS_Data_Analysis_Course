# Lesson

An RNA-Seq workflow
===================

Learning Objectives:
-------------------
### What's the goal for this lesson?

* Use a series of command line tools to execute an RNA-Seq workflow
* Automate a workflow by grouping a series of sequential commands into a script
* Modify and submit the workflow script to the cluster

## Running a Workflow

### Setting up

To get started with this lesson, ensure you are logged into the cluster and are working
in an interactive session on a compute node (single core):

```
ssh username@orchestra.med.harvard.edu
(enter password)

bsub -Is -n 6 -q interactive bash	
```

Change directories into the `unix_oct2015` directory and copy the `reference_data` folder into your project directory:

```
$ cp reference_data rnaseq_project/

```

Now move into the `rnaseq_project` directory. You should have a directory tree setup similar to that shown below. it is best practice to have all files you intend on using for your workflow present within the same directory. In our case, we have our original FASTQ files and post-trimming data generated in the previous section. We also have all reference data files that will be used in downstream analyses.

```
rnaseq_project
	├── data
	│   ├── reference_data
	│   │   └── chr1.fa
	│   │   └── chr1-hg19_genes.gtf
 	|   ├── untrimmed_fastq
	│   │   
	│   └── trimmed_fastq
	│       ├── Irrel_kd_1.qualtrim25.minlen35.fq
	│       ├── Irrel_kd_2.qualtrim25.minlen35.fq
	│       ├── Irrel_kd_3.qualtrim25.minlen35.fq
	│       ├── Mov10_oe_1.qualtrim25.minlen35.fq
	│       ├── Mov10_oe_2.qualtrim25.minlen35.fq
	│       └── Mov10_oe_3.qualtrim25.minlen35.fq
	|
	├── meta
	├── results
	└── docs
```

Without getting into the details for each step of the workflow, we first describe a general overview of the steps involved in RNA-Seq analysis:

![Workflow](../img/rnaseq-workflow.png)

1. Quality control on sequence reads
2. Trim and/or filter reads (if necessary)
3. Index the reference genome for use by STAR
4. Align reads to reference genome using STAR (splice-aware aligner)
5. Count the number of reads mapping to each gene using htseq-count
6. Statistical analysis (count normalization, linear modeling using R-based tools)


We'll first perform the commands for all the above steps (run through the workflow) for a single sample.

Next, we'll create a script for the commands and test it. 

Finally, we'll modify the script to run on the cluster.

So let's get started.

The first command is to change into our working directory:

```
     cd unix_oct2015/rnaseq-project

```

Let's load up some of the modules we need for this section: **do we need any other modules??**

```
     module load samtools
```

Create an output directory for our alignment files:

```bash
mkdir results/STAR

```

In the script, we will eventually loop over all of our files and have the cluster work on each one in parallel. For now, we're going to work on just one to set up our workflow.  To start we will use the trimmed first replicate in the Mov10 overexpression group, `Mov10_oe_1.qualtrim25.minlen35.fq` 


**NOTE: if you did not follow the last section, please execute the following command:** (this will copy over the required files into your home directory.)

```bash

cp -r /groups/hbctraining/unix_oct2015_other/trimmed_fastq data/

```

### Alignment to genome
The alignment process consists of choosing an appropriate reference genome to map our reads against, and performing the read alignment using one of several splice-aware alignment tools such as [STAR](http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635) or [TopHat2](https://ccb.jhu.edu/software/tophat/index.shtml). The choice of aligner is a personal preference and also dependent on the computational resources that are available to you.
 
For this workshop we will be using STAR (Spliced Transcripts Alignment to a Reference), an aligner designed to specifically address many of the challenges of RNAseq
data mapping, and utilizes a novel strategy for spliced alignments. STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping
speed. More details on the algorithm itself can be found in the publication linked above. Aligning reads using STAR is a two step process: 1) Create a genome index 2) Map reads to the genome.


Indexing of the reference genome has already been done for you. **You do not need to run this code**. For this step you need to provide a reference genome and an annotation file. For this workshop we are using reads that originate from a small subsection of chromosome 1 (~300,00 reads) and so we are using only chr1 as the reference genome, and have provided the appropriate indices. Depending on the size of your genome, this can take awhile. 

The basic options to **generate genome indices** using STAR as follows:


* `--runThreadN`: number of threads
* `--runMode`: genomeGenerate mode
* `--genomeDir`: /path/to/store/genome_indices
* `--genomeFastaFiles`: /path/to/FASTA_file 
* `--sjdbGTFfile`: /path/to/GTF_file
* `--sjdbOverhang`: readlength -1

```
** Do not run this**
STAR --runThreadN 5 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles chr1.fa --sjdbGTFfile chr1-hg19_genes.gtf --sjdbOverhang 99

```

The basic options for **mapping reads** to the genome using STAR is as follows:

* `--runThreadN`: number of threads
* `--readFilesIn`: /path/to/FASTQ_file
* `--genomeDir`: /path/to/genome_indices
* `--outFileNamePrefix`: prefix for all output files


More details on STAR and its functionality can be found in the [user manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), we encourage you to peruse through to get familiar with all available options.


_STAR is not available as a module on Orchestra._ To run STAR we will be using an install of the software that is available on the Orchestra cluster at `/opt/bcbio/local/bin`. Since we had previously added this location to our `$PATH` we can access the software by simply using the STAR command followed by the basic parameters described above and any additional parameters. The full command is provided below for you to copy paste into your terminal. Below, we first describe some the extra parameters we have added.

Advanced parameters:

* `--outFilterMultimapNmax`: max number of multiple alignments allowed for a read
* `--outSAMstrandField`: compatability with Cufflinks (for transcriptome assembly)
* `--outReadsUnmapped`: file format for unmapped reads
* `--outSAMtype`: output filetype (SAM default)
* `--outSAMUnmapped`: what to do with unmapped reads
* `--outSAMattributes`: specify SAM attributes in output file


```
STAR --runThreadN 6 --genomeDir /groups/hbctraining/unix_oct2015_other/reference_STAR --readFilesIn data/trimmed_fastq/Mov10_oe_1.qualtrim25.minlen35.fq  --outFileNamePrefix results/STAR/Mov10_oe_1_ --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM Unsorted --outSAMunmapped Within --outSAMattributes NH HI NM MD AS
```

#### Exercise
How many files do you see in your output directory? Using the `less` command take a look at `Mov10_oe_1_Log.final.out` and answer the following questions:  

1. How many reads are uniquely mapped?
2. How many reads map to more than 10 locations on the genome?
3. How many reads are unmapped due to read length?


#### Assess the alignment (visualization)

Index the BAM file for visualization with IGV:

    samtools index results/bam/SRR098283.trimmed.aligned.sorted.bam

**Transfer files to your laptop**

Using FileZilla, transfer the following 3 files to your local machine, 
`results/bam/SRR098283.trimmed.aligned.sorted.bam`,

`results/bam/SRR098283.trimmed.aligned.sorted.bam.bai`, 

`data/ref_genome/ecoli_rel606.fasta`

**Visualize**

* Start [IGV](https://www.broadinstitute.org/software/igv/download)
* Load the genome file into IGV using the **"Load Genomes from File..."** option under the **"Genomes"** pull-down menu.
* Load the .bam file using the **"Load from File..."** option under the **"File"** pull-down menu. *IGV requires the .bai file to be in the same location as the .bam file that is loaded into IGV, but there is no direct use for that file.*



#### Exercise I - Calling Variants from all files?

That's a lot of work, yes? But you have five more FASTQ files to go...

- Try running this workflow on a different FASTQ file. What did you have to do differently
in order to get this workflow to work?
- Remembering what commands *and* what parameters to type can be pretty daunting. What can
you do to help yourself out in this regard?
- If you were to automate this process, what additional bits of information might you need?


#### Exercise II - Automating this Workflow with a Bash Script

The easiest way for you to be able to repeat this process is to capture the steps that
you've performed in a bash script. And you've already learned how to do this in previous
lessons. So...

- Using your command history, create a script file that will repeat these commands
for you. Name your script *run_variant_call_on_file.sh*. Delete your results 
directories, and run your script. Do you get all the proper output files?

One additional command we can put in the top of the script to allow you to see what
is going on is the `set -x` bash command. This debugging tool will print out every
step before it is executed.

- Insert the debugging command in your script and re-run it. How is the output different?
If you're comfortable with how this looks and runs, then comment out this line.

- In order run this workflow on another file, you'll need to make changes. Copy this file,
giving the file a similar name, and make appropriate changes to run on another input
FASTQ file. What did you have to do differently in order to get this workflow to work?

- Knowing techniques that you've learned in previous lessons, what can we do to make this
workflow more friendly to different sets of input files?

- Again, reviewing your two scripts, are there additional commonalities across scripts
or within scripts that we could optimize?


#### Exercise III - Granting our Workflow More Flexibility

A couple of changes need to be made to make this script more friendly to both changes
in the workflow and changes in files. 

The first major change is allowing a change in the filename. Thus at the start of 
the script let's capture an input parameter that must be supplied with the script name.
This input parameter will be the name of the file we want to work on:

     fq="$1"

And we'll add a shortcut to store the location to the genome reference FASTA file:

     # location to genome reference FASTA file
     genome=data/ref_genome/ecoli_rel606.fasta

Make sure you are loading all the correct modules/tools for the script to run:
    
    # set up our software environment...
    source new-modules.sh
    module load bwa
    module load samtools
    module load bcftools

Now, walking thru the code, we can make some substitutions. To index with bwa and samtools
we can run those commands on the genome variable ($genome) so these value aren't 
static & hardcoded:

     bwa index $genome
     samtools faidx $genome

We'll keep the output paths creation, as it looks fine. (Though really, we could
put results/ in a variable and declare that at the top, so we can change where the
results will be as well. We'll leave that for an optional exercise)

     # make all of our output directories
     mkdir -p results/sai
     mkdir -p results/sam
     mkdir -p results/bam
     mkdir -p results/bcf
     mkdir -p results/vcf

In the script, it is a good idea to use echo for debugging/reporting to the screen

    echo "Processing file $fq ..."

We also need to use one special trick, to extract the base name of the file
(without the path and .fastq extension). We'll assign it
to the $base variable

    # grab base of filename for future naming
    base=$(basename $fq .fastq)
    echo "basename is $base"

Since we've already created our output directories, we can now specify all of our
output files in their proper locations. We will assign various file names to
 variables both for convenience but also to make it easier to see what 
is going on in the sommand below.

    # set up output filenames and locations
    fq=data/trimmed_fastq/$base\.fastq
    sai=results/sai/$base\_aligned.sai
    sam=results/sam/$base\_aligned.sam
    bam=results/bam/$base\_aligned.bam
    sorted_bam=results/bam/$base\_aligned_sorted.bam
    raw_bcf=results/bcf/$base\_raw.bcf
    variants=results/bcf/$base\_variants.bcf
    final_variants=results/vcf/$base\_final_variants.vcf

Our data are now staged.  We now need to change the series of command below
to use our variables so that it will run with more flexibility the steps of the 
analytical workflow>

    # Align the reads to the reference genome
    bwa aln $genome $fq > $sai

    # Convert the output to the SAM format
    bwa samse $genome $sai $fq > $sam

    # Convert the SAM file to BAM format
    samtools view -S -b $sam > $bam

    # Sort the BAM file
    samtools sort -O 'bam' -T temp.prefix $bam $sorted_bam

    # Index the BAM file for display purposes
    samtools index $sorted_bam

    # Do the first pass on variant calling by counting read coverage
    samtools mpileup -g -f $genome $sorted_bam > $raw_bcf

    # Do the SNP calling with bcftools
    bcftools call -vc -O b $raw_bcf > $variants

    # And finally, filter the SNPs for the final output
    bcftools view $variants | vcfutils.pl varFilter - > $final_variants

This new script is now ready for running:
	
	sh run_variant_call_on_file.sh <name of fastq>

#### Exercise IV - Parallelizing workflow for efficiency

To run the same script on a worker node on the cluster via the job scheduler, we need to add our **SLURM directives** at the **beginning** of the script. This is so that the scheduler knows what resources we need in order to run our job on the
compute node(s). 

Copy the `run_variant_call_on_file.sh` file and give it a new name `run_variant_call_on_file.sbatch`. Add the SLURM directives to the beginning of the file.

So the top of the file should look like:

    #!/bin/bash
    #
    #SBATCH -p serial_requeue   # Partition to submit to (comma separated)
    #SBATCH -n 1                # Number of cores
    #SBATCH -N 1                # Ensure that all cores are on one machine
    #SBATCH -t 0-1:00           # Runtime in D-HH:MM (or use minutes)
    #SBATCH --mem 100           # Memory in MB
    #SBATCH -J var_call_ecoli      # Job name
    #SBATCH -o var_call_ecoli.out       # File to which standard out will be written
    #SBATCH -e var_call_ecoli.err       # File to which standard err will be written
    #SBATCH --mail-type=ALL     # Type of email notification: BEGIN,END,FAIL,ALL
    #SBATCH --mail-user=<your-email@here.com> # Email to which notifications will be sent 

What we'd like to do is run this script on a compute node for every trimmed FASTQ -- pleasantly parallelizing our workflow. And now is where we'll use the for loop with the power of the cluster: 

    for fq in data/trimmed_fastq/*.fastq
    do
      sbatch run_variant_call_on_file.sh $fq
      sleep 1
    done

What you should see on the output of your screen would be the jobIDs that are returned
from the scheduler for each of the jobs that you submitted.

You can see their progress by using the squeue command (though there is a lag of
about 60 seconds between what is happening and what is reported).

Don't forget about the scancel command, should something go wrong and you need to
cancel your jobs.

* Change the script so that one can include an additional variable to point to
 a results directory.
