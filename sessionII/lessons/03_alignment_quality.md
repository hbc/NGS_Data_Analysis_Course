---
title: "SAM/BAM file and assessing quality "
author: "Meeta Mistry, Bob Freeman"
date: "Wednesday, October 7, 2015"
---

Approximate time:

## Learning objectives

* Evaluating the STAR aligner output files
* Understanding the standard alignment file (SAM/BAM) structure
* Using `samtools` to evaluate alignment quality 
* Visualizing alignment quality using IGV (genome browser)  


## Assessing alignment quality

After running our FASTQ files through the STAR aligner, you should have noticed a number of output files in the `~/ngs_course/rnaseq/results/STAR` directory. Let's take a quick look at some of the files that were generated and explore the content of some of them. What you should see, is that for each FASTQ file you have **5 output files** and a single tmp directory. Briefly, these files are described below:

* `Log.final.out` - a summary of mapping statistics for the sample
* `Aligned.sortedByCoord.out.bam` - the aligned reads, sorted by coordinate, in BAM format
* `Log.out` - a running log from STAR, with information about the run 
* `Log.progress.out` -  job progress with the number of processed reads, % of mapped reads etc., updated every ~1 minute
* `SJ.out.tab` - high confidence collapsed splice junctions in tab-delimited format. Only junctions supported by uniquely mapping reads are reported

## Mapping statistics

Having completed the alignment, the first thing we want to know is how well did our reads align to the reference. Rather than looking at each read alignment, it can be more useful to evaluate statistics that give a general overview for the sample. One of the output files from the STAR aligner contains mapping statistics, let's take a closer look at one of those files. We'll use the `less` command which allows us to scroll through it easily: 

	$ less Mov10_oe_1.Log.final.out
	
The log file provides information on reads that 1) mapped uniquely, 2) reads that mapped to mutliple locations and 3) reads that are unmapped. Additionally, we get details on splicing, insertion and deletion. From this file the most informative statistics include the **mapping rate and the number of multimappers**.

* As an example, a good quality sample will have **alteast 75% of the reads uniquely mapped**. Once values start to drop lower than 60% it's advisable to start troubleshooting. The lower the number of uniquely mapping reads means the higher the number of reads that are mapping to multiple locations. It is best to keep this number low because multi-mappers are not included when we start counting reads

> NOTE: The thresholds suggested above will vary depending on the organism that you are working with. Much of what is discussed here is in the context of working with human or mouse data. For example, 75% of mapped reads holds true only if the genome is good or mature. For badly assembled genomes we may not observe a high mapping rate, even if the actual sequence sample is good.


In addition to the aligner-specific summary we can also obtain quality metrics using tools like [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc). The input for RNA-SeQC can be one or more BAM files and the output consists of HTML reports and tab delimited files of metrics data. This tool can be valuable for comparing sequencing quality across different samples, but can also be run on individual samples as a means of quality control before continuing with downstream analysis. We will not be using this tool in the course, but some of the features are listed below:

* Transcript-annotated reads: Even if you have high genomic mapping rate for all samples, check to see where the reads are mapping. Ensure that there is not an unusually high number of **reads mapping to intronic regions** (~30% expected) and fewer than normally observed **mapping to exons** (~55%). A high intronic mapping suggests possible genomic DNA contamination and/or pre-mRNA. 
* Ribosomal RNA (rRNA) constitutes a large majority of the RNA species in any total RNA preparation. Despite depletion methods, you can never achieve complete rRNA removal. Even with Poly-A enrichment a small percentage of ribosomal RNA can stick to the enrichment beads non-specifically. **Excess ribosomal content (> 2%)** will normally have to be filtered out so that differences in rRNA mapped reads across samples do not affect alignment rates and skew subsequent normalization of the data. 
* GC bias and strand specificity
* Depth of coverage across transcript length


*** 

**Exercise**

Using the less command take a look at `Mov10_oe_1_Log.final.out` and answer the following questions:

1. How many reads map to more than 10 locations on the genome?
2. How many reads are unmapped due to read length?
3. What is the average mapped length per read?

***


## Alignment file format: SAM/BAM

The output we requested from the STAR aligner (using the appropriate parameters) is a BAM file. By default STAR will return a file in SAM format. BAM is a binary, compressed version of the SAM file, also known as **Sequence Alignment Map format**. The SAM file, introduced is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we will go into some features of the SAM format, the paper by [Heng Li et al](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of alignment, etc., this will change depending on the aligner being used. Each section begins with character ‘@’ followed by a two-letter record type code.  These are followed by two-letter tags and values. Example of some common sections are provided below:

```
@HD  The header line
VN: format version
SO: Sorting order of alignments

@SQ  Reference sequence dictionary
SN: reference sequence name
LN: reference sequence length
SP: species

@RG  Read group
ID: read group identifier
CN: name of sequencing center
SM: sample name

@PG  Program
PN: program name
VN: program version
```

Following the header is the **alignment section**. Each line that follows corresponds to alignment information for a single read. Each alignment line has **11 mandatory fields for essential mapping information** and a variable number of other fields for aligner specific information. 

![SAM](../img/SAM_file.png)

These fields contain information describing the read, quality of the read, and nature alignment of the read to a region of the genome. Below are example entries of alignment information for a single read (*note that in the SAM file this would be displayed tab-delimited, on a single line*): 

```
QNAME  e.g.  M00628:11:000000000-A1P5L:1:1112:26953:13136
FLAG   e.g.  163
RNAME  e.g.  CP000921
POS    e.g.  20
MAPQ   e.g.  60
CIGAR  e.g.  149M
RNEXT  e.g.  = 
PNEXT  e.g.  108
TLEN   e.g.  239
SEQ    e.g.  CCACTATGTTTTTCGATAAAAAGCTTAATAAAT
QUAL   e.g.  ?????BBBBBDBDB=?FFECFACCFFHHH>09C

```
Most of the field entries above are pretty self-explanatory, except for two which could use a bit more in-depth of an explanation. These are, the (bitwise) `FLAG` field and `CIGAR`.

### Bitwise flags explained

The `FLAG` value that is displayed can be translated into information about the mapping. The flag value corresponds to a bit set which is a sum of the binary representation of the individual flags.

There are 11 bitwise flags describing the alignment:

![bitwise](../img/bitwiseflags.png)

* For a given alignment, each of these flags are either on or off indicating the condition is true or false. 
* These flags are stored as a binary strings of length 11 instead of 11 columns of data. Value of ‘1’ indicates the flag is set.  e.g. 00100000000
* The combination of all flags are represented in the SAM specification using its hexidecimal representation
* There is no combination that can be a result of two different sums

So in our example alignment we have a bitwise flag of 163. This is the flag is actually a sum of the 4 different flags.

`163 = 1 + 2 + 32 + +128   (e.g. 00010100011)`

Which tells us that **(1)** the read is paired, **(2)** the read is mapped in a proper pair, **(32)** the read is the mate reverse strand and **(128)** the read is the second read in the pair. Bit flags are useful in storing lots of information in little space. Normally, you wouldn't be going through the SAM file manually to evaluate these numbers, but if you were curious about a particular flag and what it means, there is a [Picard utility](http://broadinstitute.github.io/picard/explain-flags.html) which can help decipher for you.

### CIGAR string explained

The CIGAR string is a sequence of letters and numbers that represent the *edits or operations* required to match the read to the reference. The letters are operations that are used to indicate which bases align to the reference (i.e. match, mismatch, deletion, insertion), and the numbers indicate the associated base lengths for each 'operation'.

![cigar](../img/cigar_strings.png)

In our example, SAM entry above the CIGAR string listed is 149M which translates to 149 matches with the reference. 

## `samtools`

[SAMtools](http://samtools.sourceforge.net/) is a tool that provides alot of functionality in dealing with SAM files. SAMtools utilities include, but are not limited to, viewing, sorting, filtering, merging, and indexing alignments in the SAM format. In this lesson we will explore a few of these utilities on our alignment files. Let's get started by loading the `samtools` module:

	module load seq/samtools/1.2

### Viewing the SAM file

Now that we have learned so much about the SAM file format, let's use `samtools` to take a quick peek at our own files. The output we had requested from STAR was a BAM file. The problem is the BAM file is binary and not human-readable. Using the `view` command within `samtools` we can easily convert the BAM into something that we can understand. You will be returned to screen the entire SAM file, and so we can either write to file, or pipe this to the `less` command so we can scroll through it.

We will do the latter (since we don't really need it for downstream analysis) and scroll through the SAM file (using the up and down arrows) to see how the fields correspond to what we expected. Adding the `-h` flag allows to also view the header.

```
$ samtools view -h results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam | less

``` 

### Summarizing and filtering the SAM file

As mentioned previously, manually reading the file line-by-line isn't very productive. It is more useful to be able to summarize across the entire file. Suppose we wanted to set a threshold on mapping quality. For example, we want to know how many reads aligned with a quality score higher than 30. To do this, we can combine the `view` command with additional flags `q 30` and `-c` (to count):

```
$ samtools view -q 30 -c results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam 

```
*How many of reads have a mapping quality of 30 or higher?*

If we wanted to only work with high quality mapped reads, we could subset these alignments and write them to a new BAM file using the `-b` flag:

```
$ samtools view -q 30 -b results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam Mov10_oe_1_Aligned_q30.bam

```

The `flagstat` command can also used to summarize a BAM file:

```
$ samtools flagstat results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam 

```

*What do you see?* The output from `flagstat` corresponds to information encoded within the 11 bitwise flags that we discussed earlier. For each flag you get the corresponding number of reads that agree. Since we are working with single end reads, many of these do not apply to us.

```
309507 + 0 in total (QC-passed reads + QC-failed reads)
9270 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
295771 + 0 mapped (95.56% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

We can apply filters to keep/remove selected reads based on where they fall within these different categories. Similar to when filtering by quality we need to use the `samtools view` command, however this time use the `-F` or `-f` flags.

These flags are combined with the hexadecimal value of the bitwise flag you are interested in. We know that 4 is the bitwise flag for unmapped reads, so to obtain the number of unmapped reads you would use `-f 4`:

```
$ samtools view -f 4 -c results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam 

```

But how do we obtain the number of reads that *are mapped*? Remember that the bitwise flags are like boolean values. If the flag exists, the condition is true. To find the number of mapped reads we need to count those reads that do not meet the condition. We do this using the capitalized F flag:

  
```
$ samtools view -F 4 -c results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam 

```

*This number should be identical to that reported on line 5 of the `flagstat` results*

### Indexing the BAM file

To perform some functions (i.e. subsetting, visualization) on the BAM file, an index is required. Indexing aims to achieve fast retrieval of alignments overlapping a specified region without going through the whole alignments. In order to index a BAM file, it must be sorted by the reference ID and then the leftmost coordinate, which can also be done with `samtools`. However, in our case we had included a parameter in our STAR alignment run so we know our BAM files are already sorted.

To index the BAM file we use the `index` command:

    $ samtools index results/STAR/Mov10_oe_1__Aligned.sortedByCoord.out.bam

This will create an index in the same directory as the BAM file, which will be identical to the input file in name but with an added extension of `.bai`.


### Subsetting the BAM file

Suppose we only wanted to look at a subset of the reads mapping to a specific location on the genome. We can extract these read alignments by using the `view` command and providing the genomic coordinates. Let's take a look at reads mapping to `chr1:200000-500000`:

	samtools view results/STAR/Mov10_oe_1__Aligned.sortedByCoord.out.bam chr1:200000-500000

*Note this will only work, if you have the index file created.* There should only be a few reads printed to screen. Use the `-c` flag to count how many entries are within the specified region.


****

**Exercise:**

1. The STAR log file for `Mov10_oe_1` indicated that there were a certain number of reads mapping to multiple locations. When this happens, one of these alignments is considered
primary and all the other alignments have the secondary alignment flag set in the SAM records. **Use `samtools` and your knowledge of [bitwise flags](https://github.com/hbc/NGS_Data_Analysis_Course/blob/master/sessionII/lessons/03_alignment_quality.md#bitwise-flags-explained) to extract the secondary reads to a file called `Mov10_oe_1_secondary_alignments.bam`.**

***


## Visualization

Another method for assessing the quality of your alignment is to visualize the alignment using a genome browser. For this course we will be using the [Integrative Genomics Viewer (IGV)(https://www.broadinstitute.org/igv/) from the Broad Institute. *You should already have this downloaded on your laptop.* IGV is an interactive tool which allows exploration of large, integrated genomic datasets. It supports a wide variety of data types, including array-based and next-generation sequence data, and genomic annotations, which facilitates invaluable comparisons.

### Transfer files

In order to visualize our alignments we will first need to move over the relevant files. We previously used FileZilla to transfer files from Orchestra to your laptop. However, there is another way to do so using the command line interface. **This option is only available for Mac and Linux users! PC users can use Filezilla.**  Similar to the `cp` command to copy there is a command that allows you to securely copy files between computers. The command is called `scp` and allows files to be copied to, from, or between different hosts. It uses ssh for data transfer and provides the same authentication and same level of security as ssh. 

First, identify the location of the _origin file_ you intend to copy, followed by the _destination_ of that file. Since the origin file is located on Orchestra, this requires you to provide remote host and login information.

The following 2 files need to be moved from Orchestra to your local machine,
 
`results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam`,

`results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam.bai` 

```
$ scp user_name@orchestra.med.harvard.edu:/home/user_name/ngs_course/rnaseq/results/Mov10_oe_1_Aligned.sortedByCoord.out.bam* /path/to/directory_on_laptop
```


### Visualize

* Start [IGV](https://www.broadinstitute.org/software/igv/download) _You should have this previously installed on your laptop_
* Load the Human genome (hg19) into IGV using the dropdown menu at the top left of your screen. _Note: there is also an option to "Load Genomes from File..." under the "Genomes" pull-down menu - this is useful when working with non-model organisms_
* Load the .bam file using the **"Load from File..."** option under the **"File"** pull-down menu. *IGV requires the .bai file to be in the same location as the .bam file that is loaded into IGV, but there is no direct use for that file.*

![IGV screenshot](../img/igv_screenshot.png)

***

**Exercise**

Now that we have done this for one sample, let's try using the same commands to perform alignment on one of the control samples. Using `Irrel_kd_1_qualtrim25.minlen35.fq` walk through the alignment commands above. Copy over the resulting BAM and index file to your laptop and upload into IGV for visualization. 

1. How does the MOV10 gene look in the control sample in comparison to the overexpression sample?
2. Take a look at a few other genes by typing into the search bar. For example, PPM1J and PTPN22. How do these genes compare? 

***

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
