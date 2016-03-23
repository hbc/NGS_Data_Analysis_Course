---
title: "Alignment with BWA"
author: "Meeta Mistry, Mary Piper"
date: "Wednesday, March 23, 2016"
---

Approximate time: 90 minutes

## Learning Objectives:

* Exploring the variant calling workflow
* Choosing appropriate BWA alignment parameters for our dataset
* Understanding alignment clean-up steps

## Variant Calling Workflow

The variant calling workflow begins with quality control and alignment, similar to the other NGS applications. Alignment is followed by alignment clean-up to prepare data for variant calling. Then, variant calling is performed, followed by filtering and annotation of the variant calls.

![var_calling_workflow](../img/variant_calling_workflow.png)

## Set-up

Login to Orchestra and start an interactive session with four cores:

```
$ bsub -Is -n 4 -q interactive bash
```

Change directories to the `ngs_course` directory:

```
$ cd ~/ngs_course
```

Create a `var-calling` directory and change directories into it:

```
$ mkdir var-calling

$ cd var-calling
```

Create folders for `data`and `results`:

```
$ mkdir data results
```

Within the `data` folder create folders for `untrimmed_fastq` and `reference_data`:

```
$ mkdir data/untrimmed_fastq data/reference_data
```

Within the `results` folder create a folder for `bwa`:

```
$ mkdir results/bwa
```

Now that we have the directory structure created, let's copy over the data to perform our quality control and alignment, including our fastq files and reference data files:

```
$ cd data

$ cp /groups/hbctraining/ngs-data-analysis2016/var-calling/raw_fastq/*fastq untrimmed_fastq/

$ cp  /groups/hbctraining/ngs-data-analysis2016/chipseq/reference_data/chr20.fa reference_data/
```


## Alignment

Choice of alignment tool is often determined by the type of NGS application being conducted. We have previously used STAR for RNA-Seq data because it is fast and optimized for aligning spliced reads. For ChIP-Seq we used Bowtie2 to align the reads because it is fast and accurate. For variant calling we will use [BWA (Burrows-Wheeler Aligner)](http://bio-bwa.sourceforge.net) for alignment. 

BWA is generally slower than Bowtie2 with similar sensitivity and both tools can perform gapped alignment for the identification of indels and can effectively map paired-end reads. However, BWA is a bit more accurate and provides information on which alignments are trustworthy. Small numbers of bad alignments can result in many false variant calls, so accuracy is paramount, and is the basis of choosing BWA.

### BWA modes

Depending on read length, BWA has different modes optimized for different sequence lengths:

- **BWA-backtrack:** designed for Illumina sequence reads up to 100bp (3-step)

- **BWA-SW:** longer sequences ranged from 70bp to 1Mbp, long-read support and split alignment

- **BWA-MEM:** share similar features to BWA-SW, but BWA-MEM is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

### Aligning reads with BWA-MEM

Change directories into the `reference_data` directory:

```
$ cd ~/ngs_course/var-calling/data/reference_data
```

Load the necessary modules for alignment and clean-up:

```
$ module load seq/samtools/1.3 seq/bwa/0.7.8  seq/picard/1.138
```

#### Creating BWA-MEM index

Similar to the other alignment tools we have used, the first step in the BWA alignment is to create an index for the reference genome. Similar to Bowtie2, BWA indexes the genome with an FM Index based on the Burrows-Wheeler Transform to keep memory requirements low for the alignment process. 

The basic options for indexing the genome using BWA are:

* `-p`: prefix for all index files

```
$ bwa index -p chr20 chr20.fa
```

#### Aligning reads with BWA-MEM

Since we have our indexes created, we can get started with read alignment. Change directories to the `bwa` folder:

```
$ cd ~/ngs_course/var-calling/data
```

We will perform alignment on our paired-end reads for sample `na12878`. Details on BWA and its functionality can be found in the [user manual](http://bio-bwa.sourceforge.net/bwa.shtml); we encourage you to peruse through to get familiar with all available options.

The basic options for aligning reads to the genome using BWA-MEM are:

* `-t`: number of threads / cores
* `-M`: mark shorter split hits as secondary (for Picard compatibility)
	

```
$ bwa mem -t 4 -M 6  \
reference_data/chr20 \   # path to genome indexes including prefix
untrimmed_fastq/na12878_1.fastq untrimmed_fastq/na12878_2.fastq \    # fastq files for paired-end reads
2> ~/ngs_course/var-calling/results/bwa/bwa.err > ~/ngs_course/var-calling/results/bwa/na12878.sam     # save standard error to file and save alignment output as `na12878.sam`

```
### Alignment clean-up

The last stage of the alignment phase is marking duplicates, and it is usually only required for variant calling. We need to find reads that are likely artifacts from the PCR amplification as they can bias variant calls.

![align_cleanup](../img/workflow_cleanup.png)

If duplicates aren't marked, then the PCR-based errors will be picked up again and again as false positive variant calls. Duplicates are easy to detect since they have the same mapping information and CIGAR string:  

![dedup1](../img/dedup_begin.png)

Marking duplicates with tools such as Picard or samblaster will result in the variant caller ignoring these PCR-based errors, and instead seeing:

![dedup1](../img/dedup_end.png)

The variant caller will be more likely to discard the error, instead of calling it as a variant.

#### Sorting SAM by coordinates
[user_manual](http://broadinstitute.github.io/picard/command-line-overview.html#Overview)

```
$ java -jar /opt/picard-1.138/SortSam.jar \
INPUT=na12878.sam \
OUTPUT=na12878_sorted.sam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=LENIENT
```

#### Marking duplicates

```
$ java -jar /opt/picard-1.138/MarkDuplicates.jar \
INPUT=na12878_sorted.sam \
OUTPUT=na12878_sorted_marked.bam \
METRICS_FILE=metrics.txt \
AS=true \
VALIDATION_STRINGENCY=LENIENT
```
#### Creating index for BAM file

```
samtools index na12878_sorted_marked.bam
```

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*




