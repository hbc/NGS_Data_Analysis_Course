Learning Objectives:
-------------------

* Use the count matrix as input to an R script for differential expression analysis
* Apply Unix commands to look at the results that are generated and extract relevant information
* Familiarize with various functional analysis tools


Differential expression analysis
-------------------

At the end of the workflow from the last lesson, our final end product was a count matrix. This is a matrix in which each row represents a gene (or feature) and each column corresponds to a sample. In our dataset, we have two sample classes (control and Mov10oe) and we want to assess the difference in expression between these groups on a gene-by-gene basis.

Intuitively, it would seem that since we know which samples belong to which group we could just compute a fold change for each gene and rank. Easy, right? Not exactly. The problem is, the gene expression that we are observing is not just a result of the differences between the groups that we are investigating, rather it is a measurement of the sum of many effects. In a given biological sample the transcriptional patterns will also be changing with respect to various extraneous factors; some that we are aware of (i.e demographic factors, batch information) and other noise that we cannot attribute to any particular source. So when we are looking for differentially expressed genes we need to try an control for these effects as best we can.


### Statistical models in R

Enter [R](https://www.r-project.org/), a software environment for statistical computing and graphics. R is widely used in the field of bioinformatics, amongst various other disciplines. It can be locally installed on almost all operating systems (and it's free!) with numerous packages available that help in increasing efficency of data handling, data manipulation and data analysis. Unfortunately, discussing the specifics about R is outside the scope of this course. However, we encourage you take a look at some of the R resources listed below if you are interested in learning more. 

R is a powerful language that can be very useful for NGS data analysis, and there are many popular tools for working with RNA-Seq count data. Some of the popular tools include [edgeR](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf), [DESeq2](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf), and [limma-voom](http://www.genomebiology.com/2014/15/2/R29). All of these tools use statistical modeling of the count data to test each gene against the null hypothesis and evaluate whether or not it is significantly differentially expressed. These methods account for within-group and between-group varaibility and are also flexible enough to allow for other coavriates you may want to account for. The details on how each tool works are described thoroughly within the vignettes.


### Running R script

In order to run R on Orchestra, let's first log on to the cluster. But this time, note the addition of `-x` in our command. A number of programs with graphical user interfaces (e.g., R, Matlab) use the X11 system which lets the program run on an Orchestra computer, but show the graphics on your desktop. To do this, you need to have an X11 server running on your desktop, and your SSH connection needs to have X11 forwarding enabled. There are different instructions provided below depending on your operating system:

**For Mac Users**
Install [Xquartz](http://xquartz.macosforge.org/landing/) and have it running on your laptop before logging on to Orchestra:

	$ ssh -x user_name@orchestra.med.harvard.edu

**For Windows Users**
Install [Xming](http://sourceforge.net/projects/xming/) and have it running on your laptop.
In PuTTY:

`Connection -> SSH -> X11 -> Enable X11 forwarding`

![puttyX11](../img/puttyssh.png)

Then login in "session" with  orchestra.med.harvard.edu.

Once you are in, start an interactive session and navigate to the `rnaseq_project` directory:

	$ bsub -Is -q interactive bash
	$ cd unix_oct2015/rnaseq_project

We will be running an R script that uses the R package [edgeR](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) to identify differentially expressed genes. This package is available through the [Bioconductor](https://www.bioconductor.org/), a repository of packages for the analysis of high-throughput genomic data. There are also a few other packages that are required to generate some additional figures.

We first need to load the R module:

	$ module load stats/R/3.2.1


You can open R by simply typing `R` at the command prompt and pressing `Enter`. You are now in the R console:

![Rconsole](../img/R_screenshot.png)


To install the packages we need we have created an R script for you to run from the command line, but normally you would do so with commands in the R console. **Since we are running a script we can exit R with:**

	q()


You should find yourself back at the shell command prompt. We will first need to copy over the installation script and setup some important _environment variables_. 

	$ cp /groups/hbctraining/unix_oct2015_other/install_libraries.R .

The next few lines will create a directory in your home folder for installing any R packages. Setting the environment variable will let R know where thos directory resides.

	$ mkdir -p ~/R/library
	$ echo 'R_LIBS_USER="~/R/library"' >  $HOME/.Renviron
	$ export R_LIBS_USER="/home/user_name/R/library"

Now you should be able to run the installation script:

	
	$ Rscript install_libraries.R
	

**The installation script may take a few minutes to run**, and you will see quite a bit of text being printed to screen. The reason for this is that R is also installing any dependencies and updating existing packages as required.If the packages installed successfully you will now be able to run the DE script. We are going to run the script from the `results` directory, so let's navigate there and create a directory for the results of our analysis. We will call the directory `diffexpression`:

	$ cd results
	$ mkdir diffexpression

First, let's copy over the script file:

	$ cp /groups/hbctraining/unix_oct2015_other/DE_script.R diffexpression/

The DE script will require as input **1) your count matrix file** and **2) a metadata file**. The count matrix we generated in the last lesson and is in the `counts` directory. The metadata file is a tab-delimited file which contains any information associated with our samples. Each row corresponds to a sample and eaach column contains some information about each sample.

> If you _didn't generate this file in class_ we have a pre-computed count matrix generated that you can use:
> 
> 	$ cp /groups/hbctraining/unix_oct2015_other/counts_STAR/Mov10_rnaseq_counts_complete.txt diffexpression



	$ cp ~/unix_oct2015/other/Mov10_rnaseq_metadata.txt diffexpression

> Once you have the files copied, take a look at the metadata using `less`.

Ok, now we're all setup to run our R script! Let's run it from within our `diffexpression` directory,

	$ cd diffexpression
	$ Rscript DE_script.R Mov10_rnaseq_counts_complete.txt Mov10_rnaseq_metadata.txt 


> How many files do you get as output from the script? There should be a few PNG files and text files. Use Filezilla or `scp` to copy the images over to your laptop and take a look what was generated. 


### Gene list exploration

There are two results files generated from `DE_script.R`, a full table and significant genes table (at FDR < 0.05). Take a look at the significant results file and see what values have been reported:

	$ head DEresults_sig_table.txt

You should have a table with 6 columns in it:

1. Gene symbols (this will not have a column name, due to the nature of the `write` function)
2. logFC
3. logCPM (average expression across all samples)
4. Log odds ratio
5. P-value
6. FDR value

Since we have the full table of values we could theoretically use that and filter the genes to our discretion. We could also increase stringency by adding in a fold change criteria. The full table is also useful for investigating genes of interest that did not appear in our significant list, and give us some insight into whether the gene missed the threshold marginally or by a landslide. 


Using `wc -l` find out how many genes are identified in the significant table? Keep in mind this is generated using the truncated dataset.

	$ wc -l DEresults_sig_table.txt

For downstream analysis, the relevant information that we will require from this results table is the gene names and the FDR value. We can cut the columns to a new file and and use that as input to some functional analaysis tools.

	cut -f1,6 DEresults_sig_table.txt > Mov10_sig_genelist.txt
  
Since the list we have is generated from analaysis on a small subset of chromsome 1, using these genes as input to downstream tools will not provide any meaningful results. As such, **we have generated a list using the full dataset for these samples and can be downloaded to your laptop via [this link](./genelist_edgeR_Mov10oe_noFC.txt).** From the full dataset analysis, 453 genes were identified as significant if they had an FDR < 0.05 _and_ a log fold change > 1.  


#### gProfiler

[gProfileR](http://biit.cs.ut.ee/gprofiler/index.cgi) is a web-based tool for the interpretation of large gene lists. The core tool takes a gene list as input and performs statistical enrichment analysis to provide interpretation to user-provided gene lists. Multiple sources of functional evidence are considered, including Gene Ontology terms, biological pathways, regulatory motifs of transcription factors and microRNAs, human disease annotations and protein-protein interactions. The user selects the organism and the sources of evidence to test. There are also additional parameters to change various thresholds and tweak the stringency to the desired level. 

![gprofiler](../img/gProfiler.png)

Take your gene list and paste it in the `Query' box. 

* Under **Options**: keep all defaults checked but for _Hierarchical Filtering_ use the pulldown to select _Best per parent_
* Choose **Show advanced options** and change the _Significance threshold_ to _Benjamini-Hochberg_
* From the functional evidence selections choose the following: Gene Ontology (biological process, molecular function), KEGG, Reactome
* Press **g:Profile!** 


> Take a look at the list of terms that appear. Do you see anything relevant, given what you knwo about this dataset? Run the analysis again but this time change the appropriate parameter to export your results to file. 


#### GeneMANIA

[GeneMANIA](http://genemania.org/) is another tool for predicting the function of your genes. Rather than looking for enrichment, the query gene set is evaluated in the context of curated functional association data and results are displayed in the form of a network. Association data include protein and genetic interactions, pathways, co-expression, co-localization and protein domain similarity. Genes are represented as the nodes of the network and edges are formed by known association evidence. The query gene set is highlighted and so you can find other genes that are related based on the toplogy in the network. This tool is more useful for smaller gene sets (< 400 genes), as you can see in the figure below our input results in a bit of a hairball that is hard to interpret.

![genemania](../img/genemania.png)

> Use the significant gene list generated from the analysis we performed in class (~40 genes) as input to GeneMANIA. Using only pathway and coexpression data as evidence, take a look at the network that results. Can you predict anything functionally from this set of genes? 


### Resources for R

* https://www.datacamp.com/courses/free-introduction-to-r
* Software Carpentry materials: http://swcarpentry.github.io/r-novice-inflammation/
* Data Carpentry materials: http://tracykteal.github.io/R-genomics/
* Materials from IQSS at Harvard: http://tutorials.iq.harvard.edu/R/Rintro/Rintro.html
* [swirl](http://swirlstats.com/): learn R interactively from within the R console
* The free "try R" class from [Code School](http://tryr.codeschool.com)
* HarvardX course ["Statistics and R for the Life Sciences"](https://courses.edx.org/courses/HarvardX/PH525.1x/1T2015/info)
