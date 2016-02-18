---
title: "Advanced plotting and data vsualization in R"
author: "Mary Piper, Meeta Mistry"
date: "Wednesday, October 28, 2015"
---

Approximate time: 60 minutes

## Learning Objectives 

* Visualizing data using basic plots in R
* Advanced plots (introducing `ggplot`)
* Exporting images to file



## Data Visualization

When we are working with large sets of numbers it can be useful to display that information graphically to develop insight. Visualization deserves an entire course of its own (there is that much to know!). During this lesson we will get you started with the basics of plotting by exploring a few features of R's base plotting package and then compare those to some more of advanced features using the popular Bioconductor package `ggplot2`.


## Basic plots in R
R has a number of built-in tools for basic graph types such as hisotgrams, scatter plots, bar charts, boxplots and much [more](http://www.statmethods.net/graphs/). Rather than going through all of different types, we will focus on the scatterplot to give you an idea of the different parameters involved in adding different features to plots when using base R functions.


### Scatterplot
A scatter plot provides a graphical view of the relationship between two sets of continuous numeric data. From our metadata file we will the `samplemeans` information and plot it against `age_in_days`, to see how mean expression changes with age. The base R function to do this is `plot(y ~ x, data)`:


	plot(samplemeans ~ age_in_days, data=metadata)


 ![scatter-1](../img/scatter-plot1.png) 

Each point represents a sample. The values on the y-axis correspond to the average expression for each sample which is dependent on the x-axis variable `age_in_days`. This plot is in it's simplest form. We can customize many features of the plot (fonts, colors, axes, titles) through [graphic options](http://www.statmethods.net/advgraphs/parameters.html).

For example, lets's start by giving our plot a title and renaming the axes. We can do that by simply adding the options `xlab`, `ylab` and `main` as arguments to the `plot()` function:

	plot(samplemeans ~ age_in_days, data=metadata, main="Expression changes with age", xlab="Age (days)", ylab="Mean expression")
	
![scatter-2](../img/scatter-plot2.png) 


We can also change the shape of the data point using `pch` and the size of the point using `cex` (to give the amount to magnify relative to the default).


	plot(samplemeans ~ age_in_days, data=metadata, main="Expression changes with age", xlab="Age (days)", ylab="Mean expression", pch="*", cex=2.0)

![scatter-3](../img/scatter-plot3.png)


We can add some color to the data points on the plot by adding `col="blue"`. Alternatively, you can sub in any another color with any of the default colors or you can [experiment with other R packages](http://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html#basic-color-specification-and-the-default-palette) to fiddle with better palettes. 

We can also add color to separate the data points by information in our data frame. For example, supppose we wanted to the data points colored by celltype. We would need to specify a vector of colours and provide the factor by which we are separating samples:

	  	plot(samplemeans ~ age_in_days, data=metadata, main="Expression changes with age", xlab="Age (days)", ylab="Mean expression", pch="*", cex=2.0, col=c("blue", "red")[celltype]])


![scatter-4](../img/scatter-plot4.png)



***

**Exercise** 


1. Change the color scheme in the scatter plot, such that it reflects the genotype of samples rather than celltype.

2. Use R help to find out how to increase the size of the text on teh axis labels.

***



## Advanced figures (`ggplot2`)


More recently, R users have moved away from base graphic options and towards a plotting package called [`ggplot2`](http://docs.ggplot2.org/) that adds a lot of functionality to the basic plots seen in the previous lesson. The syntax takes some getting used to but it's extremely powerful and flexible. We can start by re-creating some of the above plots but using ggplot functions to get a feel for the syntax.

`ggplot2` is best used on data in the `data.frame` form, so we will will work with `df` for the following figures. Let's start by loading the `ggplot2` library.

```{r}
library(ggplot2)
```

The `ggplot()` function is used to **initialize the basic graph structure**, then we add to it. The basic idea is that you specify different parts of the plot, and add them together using the `+` operator.

We will start with a blank plot and will find that you will get an error, because you need to **add layers**.

```{r, eval=FALSE}
ggplot(df) # note the error 
```

**Geometric objects** are the actual marks we put on a plot. Examples include:

* points (`geom_point`, for scatter plots, dot plots, etc)
* lines (`geom_line`, for time series, trend lines, etc)
* boxplot (`geom_boxplot`, for, well, boxplots!)

A plot **must have at least one geom**; there is no upper limit. You can add a geom to a plot using the + operator

```{r, eval=FALSE}
ggplot(df) +
  geom_point() # note what happens here
```

For a more exhaustive list on all possible geometric objects and when to use them check out [Hadley Wickham's RPubs](http://rpubs.com/hadley/ggplot2-layers). Each type of geom usually has a **required set of aesthetics** to be set, and usually accepts only a subset of all aesthetics --refer to the geom help pages to see what mappings each geom accepts. Aesthetic mappings are set with the aes() function. Examples include:

* position (i.e., on the x and y axes)
* color ("outside" color)
* fill ("inside" color) shape (of points)
* linetype
* size

To start, we will add position for the x- and y-axis since `geom_point` requires mappings for x and y, all others are optional.

```{r, fig.align='center'}
ggplot(df) +
     geom_point(aes(x = row.names(df), y= samplemeans))
```

 ![ggscatter1](../img/gg-scatter-1.png) 


Now that we have the required aesthetics, let's add some extras like color to the plot. We can `color` the points on the plot based on genotype, by specifying the column of information in our data frame:

```{r, fig.align='center'}
ggplot(df) +
  geom_point(aes(x = row.names(df), y= samplemeans, color = genotype)) 
```

 ![ggscatter1.1](../img/gg-scatter-1.1.png) 


Alternatively, we could color based on celltype. If we wanted both celltype and genotype identified on the plot, the `shape` aesthetic would allow us to assign the shape of data point. Add in `shape = celltype` to your aesthetic and see how it changes your plot:

```{r, fig.align='center'}
ggplot(df) +
  geom_point(aes(x = row.names(df), y= samplemeans, color = genotype, shape = celltype)) 
```

 ![ggscatter1.2](../img/gg-scatter-1.2.png) 


The size of the data points are quite small. We can adjust within the `geom_point()` but does not need to be included in `aes()` since we are specifying how large we want the data points, rather than mapping it to a variable.


```{r, fig.align='center'}
ggplot(df) +
  geom_point(aes(x = row.names(df), y= samplemeans, color = genotype, shape = celltype), size = rel(3.0)) 
```

 ![ggscatter1.3](../img/gg-scatter-1.3.png) 

The labels on the x-axis ticks are also quite hard to read. To do this we need to add an additional **theme layer**. The ggplot2 `theme` system handles non-data plot elements such as:

* Axis label aesthetics
* Plot background
* Facet label backround
* Legend appearance

There are built-in themes we can use, or we can adjust specific elements. For our figure we will change the x-axis labels to be plotted on a 45 degree angle with a small horizontal shift to avoid overlap.

 
```{r, fig.align='center'}
ggplot(df) +
  geom_point(aes(x = row.names(df), y= samplemeans, color = genotype, shape = celltype), size = rel(3.0)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
```

 ![ggscatter2](../img/gg-scatter-2.png) 


You may have noticed that by default ggplot has re-ordered our samples based on the alphabetical order of the row names. However, we would like to preserve the original order of our metadata dataframe within our plot. To do so we need to specify the order, by adding an additional column to our data frame:

	df$ordering <- 1:12

Now, when we plot the data we can reorder the dataframe based on the `ordering` column using the `reorder` function:

```
 ggplot(df) +
    geom_point(aes(x = reorder(row.names(df), ordering), y= samplemeans, color = genotype, shape = celltype), size = rel(3.0)) +
    theme(axis.text.x = element_text(angle=45, hjust=1))

```

 ![ggscatterordered](../img/gg-scatter-ordered.png) 



### Exporting figures to file

There are two ways in which figures and plots can be output to a file (rather than simply displaying on screen). The first (and easiest) is to export directly from the RStudio 'Plots' panel, by clicking on `Export` when the image is plotted. This will give you the option of `png` or `pdf` and selecting the directory to which you wish to save it to. 


The second option is to use R functions in the console, allowing you the flexibility to specify parameters to dictate the size and resolution of the output image.  In R’s terminology, output is directed to a particular output device and that dictates the output format that will be produced.  A device must be created or “opened” in order to receive graphical output and, for devices that create a file
on disk, the device must also be closed in order to complete the output.

Let's print our scatterplot to a pdf file format. First you need to initialize a plot using a function which specifies the graphical format you intend on creating i.e.`pdf()`, `png()`, `tiff()` etc. Within the function you will need to specify a name for your image, and the with and height (optional). This will open up the device that you wish to write to:

	pdf("figure/scatterplot.pdf")


Then we plot the image to the device, using the ggplot scatterplot that we just created. Finally, close the file using the `dev.off()` function. There are also `bmp`, `tiff`, and `jpeg` functions, though the jpeg function has proven less stable than the others.

```
 ggplot(df) +
    geom_point(aes(x = reorder(row.names(df), ordering), y= samplemeans, color = genotype, shape = celltype), size = rel(3.0)) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()
```


***

**Exercise**

1. The current axis labels default to what we gave as input to `geom_point`. We can change this by adding additional layers called `xlab()` and `ylab()` for the x- and y-axis, respectively. Add these layers to the current plot such that the x-axis is labeled "samples" and the y-axis is labeled "mean expression".
2. Use the `ggtitle` layer to add a title to your plot.
3. Export this image to a pdf file called "revised_scatterplot.pdf".

***


## Histogram

To plot a histogram we require another type of geometric object `geom_bar`, which requires a statistical transformation. Some plot types (such as scatterplots) do not require transformations, each point is plotted at x and y coordinates equal to the original value. Other plots, such as boxplots, histograms, prediction lines etc. need to be transformed, and usually has a default statistic that can be changed via the `stat_bin` argument. 

```{r, eval=FALSE}
ggplot(df) +
  geom_bar(aes(x = samplemeans))
  
```

 ![gghist1](../img/gg-hist-1.png) 


Try plotting with the default value and compare it to the plot using the binwidth values. How do they differ?

```{r, fig.align='center'}
ggplot(df) +
  geom_bar(aes(x = samplemeans), stat = "bin", binwidth=0.8)
  
```

 ![gghist2](../img/gg-hist-2.png) 



### Barplot
For a barplot, we also use the geometric object `geom_bar` except we need to change the `stat` argument to `identity` to use the actual values. Since we don't have an x variable, we need to specify the row names as our index so each sample is plotted on its own. For `fill` you can use `genotype` or `celltype` and see how the plot changes. This time we also have additional layers to specify the labels for the x- and y-axis and a main title for the plot.


```r
ggplot(df) +
  geom_bar(aes(x=row.names(df), y=samplemeans, fill=genotype), colour="black", stat="identity") +
  ggtitle('Average expression for each sample') +
  xlab('') +
  ylab('Mean expression') +
  theme(plot.title = element_text(size = rel(2.0)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)),
        axis.text.x = element_text(angle=45, vjust=0.5, hjust=0.6, size = rel(1.25)))
```

 ![ggbar-1](../figure/unnamed-chunk-22-1.png) 


## Boxplot

Now that we have all the required information for plotting with ggplot2 let's try plotting a boxplot similar to what we had done using the base plot functions at the start of this lesson. We can add some additional layers to include a plot title and change the axis labels. Explore the code below and all the different layers that we have added to understand what each layer contributes to the final graphic.


```r
ggplot(df) + 
  geom_boxplot(aes(x= genotype, y=samplemeans, fill=celltype)) + 
  ggtitle('Genotype differences in average gene expression') +
  xlab('Genotype') +
  ylab('Mean expression') +
  theme(plot.title = element_text(size = rel(2.0)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))
```

 ![ggbox-2](../figure/unnamed-chunk-21-1.png)



---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
