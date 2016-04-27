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

Before we start with variant calling, we need to set-up our directory structure, and make sure the tools are readily available. 

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

$ cp /groups/hbctraining/ngs-data-analysis2016/var-calling/raw_fastq/*fq untrimmed_fastq/

$ cp /groups/hbctraining/ngs-data-analysis2016/var-calling/reference_data/chr20.fa reference_data/
```

Now that we have the data, let's make sure that bcbio tools (`/opt/bcbio/centos/bin`) are in your PATH. First, test if you have already have them in your path:

	$ which picard
	
**If the output is `/opt/bcbio/centos/bin/picard`, then you are all set!** If you don't, then do one of the following:

**Option #1**:

	$ PATH=/opt/bcbio/centos/bin:$PATH

OR

**Option #2**, add the following line to your `.bashrc` file:

	export PATH=/opt/bcbio/centos/bin:$PATH

> *NOTE: If you would like to use the tools/programs installed outside of the bcbio set up, we have a small section at the end of this markdown which tells you how to. For today's class, please use the bcbio installations of the tools.*

## Dataset

To explore the variant calling workflow, we will be using a subset of a human WGS dataset attained from the [Genome in a Bottle Consortium (GIAB)](https://sites.stanford.edu/abms/giab). GIAB was initiated in 2011 by the National Institute of Standards and Technology "to develop the technical infrastructure (reference standards, reference methods, and reference data) to enable translation of whole human genome sequencing to clinical practice" [[1](https://sites.stanford.edu/abms/giab)].

The human WGS dataset completed by GIAB is "essentially the first complete human genome to have been extensively sequenced and re-sequenced by multiple techniques, with the results weighted and analyzed to eliminate as much variation and error as possible" [[2](https://sites.stanford.edu/abms/content/how-well-did-you-sequence-genome-nist-consortium-partners-have-answer)]". To minimize bias from any specific DNA sequencing method, the dataset was sequenced separately by 14 different sequencing experiments and 5 different platforms [[2](https://sites.stanford.edu/abms/content/how-well-did-you-sequence-genome-nist-consortium-partners-have-answer)]. The dataset acts as a 'truth set' for variation in the human genome to be used as a genotype reference set to compare variant calls against.

The source DNA, known as NA12878, was taken from a single person: the mother in a father-mother-child 'trio' (~8300 vials of DNA from a homogenized large batch of NA12878 cells for distribution to other labs). Father-mother-child 'trios' are often sequenced to utilize genetic links between family members [[2](https://sites.stanford.edu/abms/content/how-well-did-you-sequence-genome-nist-consortium-partners-have-answer)].

"The Genome in a Bottle consortium also plans to develop well-characterized whole genome reference materials from two genetically diverse groups: Asians and Ashkenazi Jews. Both reference sets will include sequenced genes from father-mother-child 'trios'" and are expected to be released in 2016 [[2](https://sites.stanford.edu/abms/content/how-well-did-you-sequence-genome-nist-consortium-partners-have-answer)].

While the sample NA12878 was sequenced at a depth of 300x, we will only be using a subset of the dataset aligning to chromosome 20. The sequencing files we will be using for NA12878 sample will have a total of ~ 4 million paired-end reads.

## QC and Alignment

In our workflow, we are going to skip over the Quality Control steps, but we will assume that we used *FastQC* to ensure there are no obvious problems with our samples and no adapter or vector contamination. Since the aligner we will use performs soft-clipping, we will skip the quality trimming step as well.

Choice of alignment tool is often determined by the type of NGS application being conducted. We have previously used STAR for RNA-Seq data because it is fast and optimized for aligning spliced reads. For ChIP-Seq we used Bowtie2 to align the reads because it is fast and accurate. For variant calling we will use [BWA (Burrows-Wheeler Aligner)](http://bio-bwa.sourceforge.net) for alignment. 

BWA is generally slower than Bowtie2 with similar sensitivity and both tools can perform gapped alignment for the identification of indels and can effectively map paired-end reads. However, BWA is a bit more accurate and provides information on which alignments are trustworthy. Small numbers of bad alignments can result in many false variant calls, so accuracy is paramount, and is the basis of choosing BWA.

### BWA modes

Depending on read length, BWA has different modes optimized for different sequence lengths:

- **BWA-backtrack:** designed for Illumina sequence reads up to 100bp (3-step)

- **BWA-SW:** designed for longer sequences ranging from 70bp to 1Mbp, long-read support and split alignment

- **BWA-MEM:** shares similar features to BWA-SW, but BWA-MEM is the latest, and is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

### Aligning reads with BWA-MEM

Change directories into the `reference_data` directory:

```
$ cd ~/ngs_course/var-calling/data/reference_data
```

#### Creating BWA-MEM index

Similar to the other alignment tools we have used, the first step in the BWA alignment is to create an index for the reference genome. Similar to Bowtie2, BWA indexes the genome with an FM Index based on the Burrows-Wheeler Transform to keep memory requirements low for the alignment process. 

The basic options for indexing the genome using BWA are:

* `-p`: prefix for all index files

```
$ bwa index -p chr20 chr20.fa
```

#### Aligning reads with BWA-MEM

Now that we have our indexes created, we can get started with read alignment. Change directories to the `data` folder:

```
$ cd ~/ngs_course/var-calling/data
```

We will perform alignment on our paired-end reads for sample `na12878`. Details on BWA and its functionality can be found in the [user manual](http://bio-bwa.sourceforge.net/bwa.shtml); we encourage you to peruse through to get familiar with all available options.

The basic options for aligning reads to the genome using BWA-MEM are:

* `-t`: number of threads / cores
* `-M`: mark shorter split hits as secondary (for Picard compatibility)

**NOTE:** BWA will soft-clip poor quality sequences from the ends of the reads by default, so we do not need to specify a parameter to perform soft clipping.

```
$ bwa mem -M -t 4  \
reference_data/chr20 \   # path to genome indexes including prefix
untrimmed_fastq/na12878_1.fq untrimmed_fastq/na12878_2.fq \    # fastq files for paired-end reads
2> ../results/bwa/bwa.err > \    # save standard error to file
../results/bwa/na12878.sam    # save alignment output to a SAM file

```
### Alignment clean-up

The last stage of the alignment phase is marking duplicates, and it is usually only required for variant calling. We need to find reads that are likely artifacts from the PCR amplification as they can bias variant calls.

![align_cleanup](../img/workflow_cleanup.png)

If duplicates aren't marked, then the PCR-based errors will be picked up again and again as false positive variant calls. Duplicates are easy to detect: since they have the same mapping information and CIGAR string:  

![dedup1](../img/dedup_begin.png)

Marking duplicates with tools such as *Picard* or *samblaster* will result in the variant caller ignoring these PCR-based errors, and instead seeing:

![dedup1](../img/dedup_end.png)

The variant caller will be more likely to discard the error, instead of calling it as a variant.

We will be using the [Picard](http://broadinstitute.github.io/picard/) suite of tools from the Broad Institute to sort the alignment SAM file and mark duplicates. The documentation for the tools and their usage and options is available in the [user_manual](http://broadinstitute.github.io/picard/command-line-overview.html#Tools).
 
#### Sorting SAM by coordinates
The *Picard* tool, `SortSam`, sorts an input SAM or BAM file by coordinate, queryname, etc. Input and output formats (SAM or BAM) are determined by the file extension.

The description of base options for the `SortSam` tool:

* `INPUT`:	The BAM or SAM file to sort. Required.
* `OUTPUT`:	The sorted BAM or SAM output file. Required.
* `SORT_ORDER`:	Sort order of output file Required. Possible values: {unsorted, queryname, coordinate, duplicate}
* `VALIDATION_STRINGENCY`: Validation stringency for all SAM files read by this program. Possible values: {STRICT, LENIENT, SILENT}
	
	**NOTE:** BWA can produce SAM records that are marked as unmapped but have non-zero MAPQ and/or non-"*" CIGAR. Typically this is because BWA found an alignment for the read that hangs off the end of the reference sequence. Picard considers such input to be invalid. In general, this error can be suppressed in Picard programs by passing VALIDATION_STRINGENCY=LENIENT or VALIDATION_STRINGENCY=SILENT [[3](https://sourceforge.net/p/picard/wiki/Main_Page/)]. 

```
$ cd ../results/bwa

$ picard SortSam \
INPUT=na12878.sam \
OUTPUT=na12878_sorted.sam \
SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=SILENT

```

#### Marking duplicates
The *Picard* tool, `MarkDuplicates`, can locate and tag duplicate reads (both PCR and optical/sequencing-driven) in a BAM or SAM file, where duplicate reads are defined as originating from the same original fragment of DNA. Explanation of the process of determining duplicate reads is provided in the [user_manual](http://broadinstitute.github.io/picard/command-line-overview.html#Tools).

The basic options for marking duplicates are:

* `INPUT`:	The BAM or SAM file to sort. Required.
* `OUTPUT`:	The sorted BAM or SAM output file. Required.
* `METRICS_FILE`: File to write duplication metrics to Required.
* `ASSUME_SORTED`: If true, assume that the input file is coordinate sorted even if the header says otherwise. Default value: false. Possible values: {true, false}
* `VALIDATION_STRINGENCY`: Validation stringency for all SAM files read by this program. Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT}

```
$ picard MarkDuplicates \
INPUT=na12878_sorted.sam \
OUTPUT=na12878_sorted_marked.bam \
METRICS_FILE=metrics.txt \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT
```
#### Creating index for BAM file

Now that we have a sorted BAM file that has duplicates marked, we would like to visualize our aligned reads in IGV. To do this, we need an index for our BAM file. As we have done in previous sessions, we will use *Samtools* to create the index:

```
samtools index na12878_sorted_marked.bam
```
---
---

## Don't want to use the bcbio installation of tools?

If you are not using the bcbio-nextgen tools you will have to load the necessary modules:

	$ module load seq/samtools/1.3 seq/bwa/0.7.8  seq/picard/1.138 
	
And, the command will be slightly different for running picard, which is a java program. Below is an example:

	$ java -jar /opt/picard-1.138/bin/picard.jar SortSam \
	INPUT=na12878.sam \
	OUTPUT=na12878_sorted.sam \
	SORT_ORDER=coordinate \
	VALIDATION_STRINGENCY=SILENT
---

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*




