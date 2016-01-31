---
title: "SAM/BAM file and assessing quality "
author: "Meeta Mistry, Bob Freeman"
date: "Wednesday, October 7, 2015"
---

Approximate time:

## Learning objectives

* Familiarizing with the standard alignment file (SAM/BAM) structure
* Using `samtools` to evaluate alignment quality 
* Visualizing alignment quality using IGV (genome browser)  

	
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

Which tells us that the read is paired, the read is mapped in a proper pair, it is the mate reverse strand and it is the second read in the pair. Bit flags are useful in storing lots of information in little space. Normally, you wouldn't be going through the SAM file manually to evaluate these numbers, but if you were curious about a particular flag and what it means, there is a [Picard utility[(http://broadinstitute.github.io/picard/explain-flags.html) which can help decipher for you.

### CIGAR string explained

The CIGAR string is a sequence of base lengths and associated ‘operations’ that are used to indicate which bases align to the reference (either a match or mismatch), are deleted, are inserted, represent introns, etc.



## `samtools`

Let's take a quick look at our alignment. To do so we first convert our BAM file into SAM format using samtools and then pipe it to the `less` command. This allows us to look at the contents without having to write it to file (since we don't need a SAM file for downstream analyses). We first need to load the samtools module:

	module load seq/samtools/1.2

```
$ samtools view -h results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam | less

```
 
Scroll through the SAM file and see how the fields correspond to what we expected.

****

**Exercise:**
On examining the SAM file

***

### Visualization

Index the BAM file for visualization with IGV:

    samtools index results/STAR/Mov10_oe_1__Aligned.sortedByCoord.out.bam

Indexing aims to achieve fast retrieval of alignments overlapping a specified region without going through the whole alignments. BAM must be sorted by the reference ID and then the leftmost coordinate before indexing

**Transfer files to your laptop using the command line**

We previously used FileZilla to transfer files from Orchestra to your laptop. However, there is another way to do so using the command line interface. _This option is only available for Mac and Linux users! PC users can use Filezilla._  Similar to the `cp` command to copy there is a command that allows you to securely copy files between computers. The command is called `scp` and allows files to be copied to, from, or between different hosts. It uses ssh for data transfer and provides the same authentication and same level of security as ssh. 

First, identify the location of the _origin file_ you intend to copy, followed by the _destination_ of that file. Since the origin file is located on Orchestra, this requires you to provide remote host and login information.

The following 2 files need to be moved from Orchestra to your local machine,
 
`results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam`,

`results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam.bai` 

```
$ scp user_name@orchestra.med.harvard.edu:/home/user_name/unix_oct2015/rnaseq_project/results/Mov10_oe_1_Aligned.sortedByCoord.out.bam* /path/to/directory_on_laptop
```


**Visualize**

* Start [IGV](https://www.broadinstitute.org/software/igv/download) _You should have this previously installed on your laptop_
* Load the Human genome (hg19) into IGV using the dropdown menu at the top left of your screen. _Note: there is also an option to "Load Genomes from File..." under the "Genomes" pull-down menu - this is useful when working with non-model organisms_
* Load the .bam file using the **"Load from File..."** option under the **"File"** pull-down menu. *IGV requires the .bai file to be in the same location as the .bam file that is loaded into IGV, but there is no direct use for that file.*

![IGV screenshot](../img/igv_screenshot.png)

#### Exercise
Now that we have done this for one sample, let's try using the same commands to perform alignment on one of the control samples. Using `Irrel_kd_1_qualtrim25.minlen35.fq` walk through the alignment commands above. Copy over the resulting BAM and index file to your laptop and upload into IGV for visualization. 

1. How does the MOV10 gene look in the control sample in comparison to the overexpression sample?
2. Take a look at a few other genes by typing into the search bar. For example, PPM1J and PTPN22. How do these genes compare? 


---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
