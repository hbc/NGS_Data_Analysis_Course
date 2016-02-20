Learning Objectives:
-------------------

*  determine how functions are attributed to genes using Gene Ontology terms
*  understand the theory of how functional enrichment tools yield statistically enriched functions or interactions
*  explore functional analysis tools
*  be able to discuss functional class scoring and co-expresssion clustering

# Functional analysis using gene lists

The output of RNA-Seq differential expression analysis is a list of significant differentially expressed genes (DEGs). These DEGs can be explored in greater detail to learn more about the roles these genes may be playing in the condition of interest, by performing any of the following analyses:

- determining whether there is enrichment of known biological functions or interactions, pathways, or networks
- identifying genes of novel pathways or networks by grouping genes together based on similar trends
- understanding global changes in gene expression by visualizing all genes being significantly up- or down-regulated.

Generally for any differential expression analysis, it is useful to investigate functional enrichment and pathways associated with the DEGs using freely available web-based tools.  While tools for functional analysis use a variety of techniques, there are three main types: over-representation analysis, functional class scoring, and pathway topology [[4](../../resources/pathway_tools.pdf)].

![Pathway analysis tools](../img/pathway_analysis.png)

## Over-representation analysis
There are a plethora of functional enrichment tools available to choose from; however, many of these tools query databases with information about gene function and interactions. The ability of these tools to query databases for gene function is due to the use of a consistent vocabulary by independent databases to describe gene function. This vocabulary was established by the Gene Ontology project, and the words in the vocabulary are referred to as Gene Ontology (GO) terms. 

### Gene Ontology project

"The Gene Ontology project is a collaborative effort to address the need for consistent descriptions of gene products across databases" [[1](geneontology.org/page/documentation)]. The [Gene Ontology Consortium](http://geneontology.org/page/go-consortium-contributors-list) maintains the GO terms, and these GO terms are incorporated into gene annotations in many of the popular repositories for animal, plant, and microbial genomes. Tools that investigate **enrichment of biological functions or interactions** can query these databases for GO terms associated with a list of genes to determine whether any GO terms associated with particular functions or interactions are enriched in the gene set. Therefore, to best use and interpret the results from these functional analysis tools, it is helpful to have a good understanding of the GO terms themselves.

### GO terms

#### GO Ontologies

To describe the roles of genes and gene products, GO terms are organized into three independent controlled vocabularies (ontologies) in a species-independent manner: 

- **Biological process:** refers to the biological role involving the gene or gene product, and could include "transcription", "signal transduction", and "apoptosis". A biological process generally involves a chemical or physical change of the starting material or input.
- **Molecular function:** represents the biochemical activity of the gene product, such activities could include "ligand", "GTPase", and "transporter". 
- **Cellular component:** refers to the location in the cell of the gene product. Cellular components could include "nucleus", "lysosome", and "plasma membrane".

"The relationships between a gene product to biological process, molecular function and cellular component are one-to-many, reflecting the biological reality that a particular protein may function in several processes, contain domains that carry out diverse molecular functions, and participate in multiple alternative interactions with other proteins, organelles or locations in the cell" [[2](go.pdf)]. Therefore, a single gene product can be associated with many GO terms. Each GO term has a term name (e.g. DNA repair) and a unique term accession number (GO:0005125).

#### GO term hierarchy

Some gene products are well-researched, with vast quantities of data available regarding their biological processes and functions. However, other gene products have very little data available about their roles in the cell. 

For example, the protein, "p53", would contain a wealth of information on it's roles in the cell, whereas another protein might only be known as a "membrane-bound protein" with no other information available. 

The GO ontologies were developed to describe and query biological knowledge with differing levels of information available. To do this, GO ontologies are loosely hierarchical, ranging from general, 'parent', terms to more specific, 'child' terms. The GO ontologies are "loosely" hierarchical since 'child' terms can have multiple 'parent' terms.

Some genes with less information may only be associated with general 'parent' terms or no terms at all, while other genes with a lot of information have many terms.

![Nature Reviews Cancer 7, 23-34 (January 2007)](../img/go_heirarchy.jpg)

[Tips for working with GO terms](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003343)

### Hypergeometric testing

In a set of genes, the frequency of GO terms can be determined, and the comparison of frequencies between a gene list & a “background” set will inform us about the over- or under- representation of the GO terms. This type of testing can inform us about over- or under-representation of other entities such as particular motifs or pathways too.

![go_frequencies](../img/go_freq.png)

To determine whether GO terms (or motifs and pathways) are over- or under-represented, you can determine the probability of having a certain number of genes associated with specific GO terms for the size of the gene list based on the background set. The background dataset can be all genes in genome for your organism or you can select your own background to use.

For example, let's suppose there are 13,000 total genes in the honeybee genome and 85 genes are associated with the GO term "DNA repair". In your gene list, there are 50 genes associated with "DNA repair" out of 1,000 genes in gene list. 

By comparing the ratios, 85/13,000 in "background" dataset and 50/1,000 in your gene list, it's evident that the GO term "DNA repair" is over-represented in your dataset.

To determine whether a GO term or pathway is significantly over- or under-represented, tools often perform **hypergeometric testing**. "The hypergeometric distribution is a discrete probability distribution that describes the probability of k successes in n draws, without replacement, from a finite population of size N that contains exactly K successes, wherein each draw is either a success or a failure" [[3](https://en.wikipedia.org/wiki/Hypergeometric_distribution)]. 

Therefore, using our example, the hypergeometric distribution describes the probability of 50 genes (k) being associated with "DNA repair", for all genes in our gene list (n=1,000), from a population of all of the genes in entire genome (N=13,000) which contains 85 genes (K) associated with "DNA repair".

The calculation of probability of k successes follows the formula:

![hypergeo](../img/hypergeo.png) 


### gProfiler

[gProfileR](http://biit.cs.ut.ee/gprofiler/index.cgi) is a web-based tool for the interpretation of large gene lists. The core tool takes a gene list as input and performs statistical enrichment analysis using hypergeometric testing to provide interpretation to user-provided gene lists. Multiple sources of functional evidence are considered, including Gene Ontology terms, biological pathways, regulatory motifs of transcription factors and microRNAs, human disease annotations and protein-protein interactions. The user selects the organism and the sources of evidence to test. There are also additional parameters to change various thresholds and tweak the stringency to the desired level. 

![gprofiler](../img/gProfiler.png)

You can use gProfiler for a wide selection of organisms, and the tool accepts your gene list as input. If your gene list is ordered (e.g. by padj. values), then gProfiler will take the order of the genes into account when outputting enriched terms or pathways.

In addition, a large number (70%) of the functional annotations of GO terms are determined using _in silico_ methods to infer function from electronic annotation (IEA). While these annotations can offer valuable information, the information is of lower confidence than experimental and computational studies, and these functional annotations can be easily filtered out. 

The color codes in the gProfiler output represent the quality of the evidence for the functional annotation. For example, weaker evidence is depicted in blue, while strong evidence generated by direct experiment is shown with red or orange. Similar coloring is used for pathway information, with well-researched pathway information shown in black, opposed to lighter colors. Grey coloring suggests an unknown gene product or annotation [gProfiler_paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1933153/).

Also, due to the hierarchical structure of GO terms, you may return many terms that seem redundant since they are child and parent terms. gProfiler allows for hierarchical filtering, returning only the best term per parent term.


Take your ordered gene list and paste it in the `Query' box. 

* Under **Options**: keep all defaults checked but check _ordered_query_ and for _Hierarchical Filtering_ use the pulldown to select _Best per parent_
* From the functional evidence selections choose the following: Gene Ontology (biological process, molecular function), [KEGG](http://nar.oxfordjournals.org/content/44/D1/D457.full.pdf), and [Reactome](http://www.reactome.org).
* Press **g:Profile!** 


> Take a look at the list of terms that appear. Do you see anything relevant, given what you know about this dataset? Run the analysis again but this time change the appropriate parameter to export your results to file. 

#### gProfiler in R

While the web interface for gProfiler is a bit more intuitive to understand, we don't actually need to leave R to run gProfiler:

```
### Functional analysis of MOV10 Knockdown using gProfileR (some of these are defaults; check help pages) 


library(gProfileR)

gprofiler_results_kd <- gprofiler(query = sigKD, 
                               organism = "hsapiens",
                               ordered_query = T, 
                               exclude_iea = F, 
                               max_p_value = 0.05, 
                               max_set_size = 0,
                               correction_method = "fdr",
                               hier_filtering = "none", 
                               domain_size = "annotated",
                               custom_bg = "")
                               
### Functional analysis of MOV10 Overexpression using gProfileR 

gprofiler_results_oe <- gprofiler(query = sigOE, 
                               organism = "hsapiens",
                               ordered_query = T, 
                               exclude_iea = F, 
                               max_p_value = 0.05, 
                               max_set_size = 0,
                               correction_method = "fdr",
                               hier_filtering = "none", 
                               domain_size = "annotated",
                               custom_bg = “")


```

Let's save the gProfiler results to file:

```
## Write results to file

write.table(gprofiler_results_kd, 
            file.path(resultsDir, 'gprofiler_MOV10_kd.txt'),          
            sep="\t", quote=F, row.names=F)
            
write.table(gprofiler_results_oe, 
            file.path(resultsDir, 'gprofiler_MOV10_oe.txt'),          
            sep="\t", quote=F, row.names=F)
```

Now, extract only the lines in the gProfiler results with GO term accession numbers for downstream analyses:

```
## Extract GO IDs for downstream analysis

allterms_kd <- gprofiler_results_kd$term.id

GOs_kd <- allterms[grep('GO:', allterms_kd)]

allterms_oe <- gprofiler_results_oe$term.id

GOs_oe <- allterms[grep('GO:', allterms_oe)]
```

### REVIGO

[REVIGO](http://revigo.irb.hr/) is a web-based tool that can take our list of GO terms, collapse redundant terms by emantic similarity, and summarize them graphically. 


![REVIGO_input](../img/revigo_input.png)

Copy and paste the GO ids from GOs_oe into the search box, and submit.

![REVIGO_output](../img/revigo_output.png)


## Functional class scoring tools
Functional class scoring tools, such as [GSEA](http://software.broadinstitute.org/gsea/index.jsp), use the gene-level statistics from the differential expression results for pathway analysis. These gene-level statistics for all genes in a given pathway are used to generate a single pathway-level statistic. The statistical significance of the pathway statistic is determined either by only considering genes in the given pathway or by comparing genes in a given pathway to genes not in the pathway.

## Pathway topology tools
Pathway topology-based methods utilize the number and type of interactions between gene product (our DE genes) and other gene products to identify pathways associated with the condition of interest.

### GeneMANIA

[GeneMANIA](http://genemania.org/) is another tool for predicting the function of your genes. Rather than looking for enrichment, the query gene set is evaluated in the context of curated functional association data and results are displayed in the form of a network. Association data include protein and genetic interactions, pathways, co-expression, co-localization and protein domain similarity. Genes are represented as the nodes of the network and edges are formed by known association evidence. The query gene set is highlighted and so you can find other genes that are related based on the toplogy in the network. This tool is more useful for smaller gene sets (< 400 genes), as you can see in the figure below our input results in a bit of a hairball that is hard to interpret.

![genemania](../img/genemania.png)

> Use the significant gene list generated from the analysis we performed in class as input to GeneMANIA. Using only pathway and coexpression data as evidence, take a look at the network that results. Can you predict anything functionally from this set of genes? 

### Co-expression clustering

Co-expression clustering is often used to identify genes of novel pathways or networks by grouping genes together based on similar trends. Co-clustering tools are useful in identifying genes in a pathway, when their participation in a pathway and/or the pathway itself is unknown. These tools cluster genes with similar expression patterns across conditions or in a time-course experiment, creating a "refined" gene list.

You can visualize co-expression clustering using heatmaps, which should be viewed as suggestive only; serious classification of genes needs better methods.  

The way the tools perform clustering is by taking the entire expression matrix and computing pair-wise co-expression values. A network is then generated from which we explore the topology to make inferences on gene co-regulation.

The [WGCNA](http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork ) package (in R) provides a more sophisticated method for co-expression clustering.

## Resources for functional analysis

* g:Profiler - http://biit.cs.ut.ee/gprofiler/index.cgi 
* DAVID - http://david.abcc.ncifcrf.gov/tools.jsp 
* GeneMANIA - http://www.genemania.org/
* GenePattern -  http://www.broadinstitute.org/cancer/software/genepattern/ (need to register)
* WebGestalt - http://bioinfo.vanderbilt.edu/webgestalt/ (need to register)
* AmiGO - http://amigo.geneontology.org/amigo
* ReviGO (visualizing GO analysis, input is GO terms) - http://revigo.irb.hr/ 
* WGCNA - http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork
* GSEA - http://software.broadinstitute.org/gsea/index.jsp
