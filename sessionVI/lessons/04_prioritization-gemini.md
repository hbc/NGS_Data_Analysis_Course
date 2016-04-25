---
title: "Variant priotization with Gemini"
author: "Meeta Mistry, Mary Piper"
date: "Friday, April 29th, 2016"
---

Approximate time: 90 minutes

## Learning Objectives:

* Using annotations information to filter out important variants
* Learning how to use GEMINI, a framework for exploring variant information



## Prioritizing variants 

Now we have annotations for all of our variants, but how do we easily sift through and find the important ones? To tackle this problem we look to tools beyond your text editor (simple shell scripts) and excel. Aaron Quinlan’s lab at University of Utah has been developing a framework called [GEMINI (GEnome MINIng)](https://github.com/arq5x/gemini) for quite some time now. 

GEMINI is a tool that helps turn those giant, sparse VCF variant matrices (millions of rows, thousands of columns) into a simple, accessible database. Within the database GEMINI annotates with just about everything out there. ENCODE, OMIM, dbSNP… *plus* internal annotations like regions of interest, candidate genes, etc. The resulting framework supports an **interactive exploration of variant information** in the context of numerous third-party genomic annotations.


<img src="../img/Gemini.png" width="600">


Let's start by loading our VCF file into the database. This command assumes that the VCF has been pre-annotated with snpEff as pecified with `-t`. While loading the database, GEMINI computes many additional population genetics statistics that support downstream analyses.

	$ gemini load -v results/annotation/na12878_q20_annot_snpEff.vcf -t snpEff \
       results/annotation/na12878_q20.db

### Assembling a query

To explore variants GEMINI, we need to use SQL (Structured Query Language) to create simple, powerful queries based on annotations, genotypes or a combination of both. It will take some time to get used to the language but once you have the hang of it, you‘ll see how powerful it is.

Below is an example query used to demonstrate the structure and what each of the components represent. To begin you will need to use the `gemini query` command to **ask GEMINI to report to us data from the database that matches the criteria we provide**.

<img src="../img/gemini-query1.png" width="400">

The query itself is written as SQL statements using `select`, `from` and `where`. Remember that we are querying a database, and databases most often contains one or more **tables**. Each table is identified by a name, and contain some number of **fields/columns** and **records (rows)** with data. 

When we `select` we are identifying the **fields** we are interested in retrieving. In our case this is the chromsome name, start and end coordinates. 

<img src="../img/gemini-query2.png" width="400">

The `from` is used to identify the **name of the table** we wish to query from. Here, we are querying from a table called `variants`.

<img src="../img/gemini-query3.png" width="400">

The `variants` table is one of many available in GEMINI. It contains a wide array of information on variants including the snpEff annotations, disease phenotype information, ENCODE data and much, much more. To find out more about the different tables and types of information available in GEMINI, take a look at the documentation pertaining to [database schema](http://gemini.readthedocs.org/en/latest/content/database_schema.html).

 <img src="../img/gemini-table.png" width="700"> 
 
Finally, there is also a `where` clause which is used to **filter records/rows**. We can add criteria, that each row must satisfy in order for it to be retrieved. In our case the rows are variants, and in our example we are querying for variants in which the `type` field is equal to `snp.`

<img src="../img/gemini-query4.png" width="400">

The final touches to the command involve wrapping the entire statement in double quoutation marks and then adding the name of the database that we wish to query.


 
***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*




