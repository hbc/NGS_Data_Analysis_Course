---
title: "The Shell: Searching and Redirection"
author: "Sheldon  McKay, Bob Freeman, Mary Piper"
date: "January 22, 2016"
---

Approximate time: 60 minutes

## Learning Objectives

* learn how to search for characters or patterns in a text file using the `grep` command
* learn about output redirection
* learn how to use the pipe (`|`) character to chain together commands

## Searching files

We went over how to search within a file using `less`. We can also
search within files without even opening them, using `grep`. Grep is a command-line
utility for searching plain-text data sets for lines matching a string or regular expression (regex).
Let's give it a try!

Suppose we want to see how many reads in our file `Mov10_oe_1.subset.fq` are really bad, with 10 consecutive Ns  
Let's search for the string NNNNNNNNNN: 

`$ cd ~/ngs_course/unix_lesson/raw_fastq`

`$ grep NNNNNNNNNN Mov10_oe_1.subset.fq`

We get back a lot of lines.  What if we want to see the whole fastq record for each of these reads?
We can use the '-B' and '-A' arguments for grep to return the matched line plus one before (-B 1) and two
lines after (-A 2). Since each record is four lines and the second line is the sequence, this should
return the whole record.

`$ grep -B1 -A2 NNNNNNNNNN Mov10_oe_1.subset.fq`

for example:

	@HWI-ST330:304:H045HADXX:1:1101:1111:61397
	CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
	+
	@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
	--
	@HWI-ST330:304:H045HADXX:1:1101:1106:89824
	CACAAATCGGCTCAGGAGGCTTGTAGAAAAGCTCAGCTTGACANNNNNNNNNNNNNNNNNGNGNACGAAACNNNNGNNNNNNNNNNNNNNNNNNNGTTGG
	+
	?@@DDDDDB1@?:E?;3A:1?9?E9?<?DGCDGBBDBF@;8DF#########################################################

****
**Exercise**

1. Search for the sequence CTCAATGA in `Mov10_oe_1.subset.fq`.
In addition to finding the sequence, have your search also return
the name of the sequence.
2. Search for that sequence in all Mov10 replicate fastq files.

****

## Redirection

We're excited we have all these sequences that we care about that we
just got from the FASTQ files. That is a really important motif
that is going to help us answer our important question. But all those
sequences just went whizzing by with grep. How can we capture them?

We can do that with something called "redirection". The idea is that
we're redirecting the output from the terminal (all the stuff that went
whizzing by) to something else. In this case, we want to print it
to a file, so that we can look at it later.

The redirection command for putting something in a file is `>`

Let's try it out and put all the sequences that contain 'NNNNNNNNNN'
from all the files in to another file called `bad_reads.txt`.

`$ grep -B1 -A2 NNNNNNNNNN Mov10_oe_1.subset.fq > bad_reads.txt`

The prompt should sit there a little bit, and then it should look like nothing
happened. But type `ls`. You should have a new file called `bad_reads.txt`. Take
a look at it and see if it has what you think it should.

If we use '>>', it will append to rather than overwrite a file.  This can be useful for saving more than one search, for example:
    
`$ grep -B1 -A2 NNNNNNNNNN Mov10_oe_2.subset.fq >> bad_reads.txt`

Since our `bad_reads.txt` file isn't a raw_fastq file, we should move it to a different location within our directory. We decide to create a new folder called `other`, and move the `bad_reads.txt` to this `other` folder using the command `mv`. 

`$ mkdir ../other/`

`$ mv bad_reads.txt ../other/`

There's one more useful redirection command that we're going to show, and that's
called the pipe command, and it is `|`. It's probably not a key on
used very much on your keyboard. What `|` does is take the output that
scrolling by on the terminal and then can run it through another command.
When it was all whizzing by before, we wished we could just slow it down and
look at it, like we can with `less`. Well it turns out that we can! We pipe
the `grep` command to `less`

`$ grep -B1 -A2 NNNNNNNNNN Mov10_oe_1.subset.fq | less`

Now we can use the arrows to scroll up and down and use `q` to get out.

We can also do something tricky and use the command `wc`. `wc` stands for
*word count*. It counts the number of lines or characters. So, we can use
it to count the number of lines we're getting back from our `grep` command.
And that will magically tell us how many sequences we're finding.

`$ grep NNNNNNNNNN Mov10_oe_1.subset.fq | wc`

This command tells us the number of lines, words and characters in the file. If we
just want the number of lines, we can use the `-l` flag for `lines`.

`$ grep NNNNNNNNNN Mov10_oe_1.subset.fq | wc -l`

Redirecting is not super intuitive, but it's powerful for stringing
together these different commands, so you can do whatever you need to do.

The philosophy behind these commands is that none of them
really do anything all that impressive. BUT when you start chaining
them together, you can do some really powerful things 
efficiently. If you want to be proficient at using the shell, you must
learn to become proficient with the pipe and redirection operators:
`|`, `>`, `>>`.

##Practice with searching and redirection

Finally, let's use the new tools in our kit and a few new ones to examine our gene annotation file, **chr1-hg19_genes.gtf**, which we will be using later to find the genomic coordinates of all known exons on chromosome 1.

`$ cd ~/ngs_course/unix_lesson/reference_data/`

Let's explore our `chr1-hg19_genes.gtf` file a bit. What information does it contain?

`$ less chr1-hg19_genes.gtf`

	chr1    unknown exon    14362   14829   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
	chr1    unknown exon    14970   15038   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
	chr1    unknown exon    15796   15947   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
	chr1    unknown exon    16607   16765   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
	chr1    unknown exon    16858   17055   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";

The columns in the gtf file contain the genomic coordinates of gene features (exon, start_codon, stop_codon, CDS) and the gene_names, transcript_ids and protein_ids (p_id) associated with these features. Note that sometimes an exon can be associated with multiple different genes and/or transcripts. For example, 

`$ grep FAM138* chr1-hg19_genes.gtf | head -n 5`

This search returns two different genes, FAM138A and FAM138F, that contain the same exon.

`$ grep PLEKHN1 chr1-hg19_genes.gtf | head -n 5`

This search returns two different transcripts of the same gene, NM_001160184 and NM_032129, that contain the same exon.

Now that we know what type of information is inside of our gtf file, let's explore our commands to answer a simple question about our data. Let's find how many total exons are present on chromosome 1 using our **gtf** file, `chr1-hg19_genes.gtf`. 

To determine the number of total exons on chromosome 1, we are going to perform a series of steps:
	
	1. Subset the dataset to only include the feature type and genomic location information
	2. Extract only the genomic coordinates of exon features
	3. Remove duplicate exons
	4. Count the total number of exons
	
####Subsetting dataset
Since we only need the feature type and the genomic location information to find the total number of exons, we should only keep columns 1, 3, 4, 5, and 7. 

`cut` is a program that will extract columns from files.  It is a very good command to know.  Let's first try out the `cut` command on a small dataset (just the first 5 lines of chr1-hg19_genes.gtf) to make sure we have the command correct:

`$ head -n 5 chr1-hg19_genes.gtf | cut -f1,3,4,5,7`
   
`-f1,3,4,5,7` means to cut these fields (columns) from the dataset.  

	chr1	exon	14362	14829	-
	chr1	exon	14970	15038	-
	chr1	exon	15796	15947	-
	chr1	exon	16607	16765	-
	chr1	exon	16858	17055	-

The `cut` command assumes our data columns are separated by tabs (i.e. tab-delimited). The `chr1-hg19_genes.gtf` is a tab-delimited file, so the default `cut` command works for us. However, data can be separated by other types of delimiters. Another common delimiter is the comma, which separates data in comma-separated value (csv) files. If your data is not tab delimited, there is a `cut` command argument (-d) to specify the delimiter.

Our output looks good, so let's cut these columns from the whole dataset (not just the first 5 lines) and save it as a file, '**chr1-hg19genes_cut**':

`$ cut -f1,3,4,5,7 chr1-hg19_genes.gtf > chr1-hg19genes_cut`

Check the cut file to make sure that it looks good using `less`. 

####Extracting genomic coordinates of exon features
We only want the exons (not CDS or start_codon features), so let's use `grep` to only keep the exon lines and save to file, '**chr1_exons**:

`$ grep exon chr1-hg19genes_cut > chr1_exons`

#### Removing duplicate exons
Now, we need to remove those exons that show up multiple times for different genes or transcripts.    

We can use some new tools `sort` and `uniq` to extract only those unique exons. `uniq` is a command that will omit repeated adjacent lines of data if they are exactly the same. Therefore, we need to sort our data by genomic coordinates first to make sure that all matching exons are adjacent to each other. 

We can use the `sort` command with the `-k` option for sort to specify which column(s) to sort on.  Note that this does something similar to cut's `-f`.

`$ sort -k3,4 chr1_exons | uniq`

####Counting the total number of exons
Now, to count how many unique exons are on chromosome 1, we need to pipe the output to `wc -l`:

`$ sort -k3,4 chr1_exons | uniq | wc -l`
    

****

**Exercise**

1. How could have you have determined the number of total exons by combining all of the previous commands (starting with the original chr1-hg19_genes.gtf), into a single command (no intermediate files) using pipes?

2. There is an argument for the `sort` command that will only keep unique lines of data. Determine the number of total exons without using the `uniq` command.

3. There is an argument for the `uniq` command that will count the number of occurrences of non-unique exons. Use the `uniq` command to count the number of duplicated exons and determine the most occurrences of an exon in the dataset.
****


## Where can I learn more about the shell?

- Software Carpentry tutorial - [The Unix shell](http://software-carpentry.org/v4/shell/index.html)
- The shell handout - [Command Reference](http://files.fosswire.com/2007/08/fwunixref.pdf)
- [explainshell.com](http://explainshell.com)
- http://tldp.org/HOWTO/Bash-Prog-Intro-HOWTO.html
- man bash
- Google - if you don't know how to do something, try Googling it. Other people
have probably had the same question.
- Learn by doing. There's no real other way to learn this than by trying it
out.  Write your next paper in nano (really emacs or vi), open pdfs from
the command line, automate something you don't really need to automate.


**Commands, options, and keystrokes covered in this lesson**

```bash
grep
> (output redirection)
>> (output redirection, append)
| (output redirection, pipe)
wc
cut
sort
uniq
```

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*

* *Adapted from the lesson by Tracy Teal. Contributors: Paul Wilson, Milad Fatenejad, Sasha Wood, and Radhika Khetani for Software Carpentry (http://software-carpentry.org/)*

