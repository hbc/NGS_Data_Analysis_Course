---
title: "Advanced plotting and data vsualization in R"
author: "Mary Piper, Meeta Mistry"
date: "Wednesday, October 28, 2015"
---

Approximate time: 90 minutes

## Learning Objectives 

* Visualizing data using basic plots in R
* Advanced plots (introducing `ggplot`)
* Exporting images to file



## Data Visualization

When we are working with large sets of numbers it can be useful to display that information graphically to gain more insight. Visualization deserves an entire course of its own (there is that much to know!). During this lesson we will get you started with the basics of plotting by exploring a few features of R's base plotting package and then compare those to some more of advanced features using the popular Bioconductor package `ggplot2`.


## Basic plots in R
R has a number of built-in tools for basic graph types such as hisotgrams, scatter plots, bar charts, boxplots and much [more](http://www.statmethods.net/graphs/). Rather than going through all of different types, we will focus on the scatterplot to give you an idea of the different parameters involved in changing features to R base plots.

### Scatterplot
A scatter plot provides a graphical view of the relationship between two sets of continuous numeric data. From our new_metadata file we will take the `samplemeans` column and plot it against `age_in_days`, to see how mean expression changes with age. The base R function to do this is `plot(y ~ x, data)`:


	plot(samplemeans ~ age_in_days, data=new_metadata)


 ![scatter-1](../img/scatter-plot1.png) 

Each point represents a sample. The values on the y-axis correspond to the average expression for each sample which is dependent on the x-axis variable `age_in_days`. This plot is in its simplest form, we can customize many features of the plot (fonts, colors, axes, titles) through [graphic options](http://www.statmethods.net/advgraphs/parameters.html).

For example, let's start by giving our plot a title and renaming the axes. We can do that by simply adding the options `xlab`, `ylab` and `main` as arguments to the `plot()` function:

```
plot(samplemeans ~ age_in_days, data=new_metadata, main="Expression changes with age", xlab="Age (days)", 
	ylab="Mean expression")
```	
	
![scatter-2](../img/scatter-plot2.png) 


We can also change the **shape of the data point using the `pch`** option and the **size of the data points using `cex`** (specifying the amount to magnify relative to the default).

```
plot(samplemeans ~ age_in_days, data=new_metadata, main="Expression changes with age", xlab="Age (days)", 
	ylab="Mean expression", pch="*", cex=2.0)
```

![scatter-3](../img/scatter-plot3.png)


We can also add some **color to the data points** on the plot by adding `col="blue"`. Alternatively, you can sub in any of the default colors or you can [experiment with other R packages](http://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html#basic-color-specification-and-the-default-palette) to fiddle with better palettes. 

We can also add color to **separate the data points by information** in our data frame. For example, supppose we wanted to the data points colored by celltype. We would need to specify a vector of colours and provide the factor by which we are separating samples:

```
plot(samplemeans ~ age_in_days, data=new_metadata, main="Expression changes with age", xlab="Age (days)", 
	ylab="Mean expression", pch="*", cex=2.0, col=c("blue", "green")[celltype]])
```

![scatter-4](../img/scatter-plot4.png)

The last thing this plot needs is a legend describing the color scheme. It would be great if it created one for you by default, but with R base functions unfortunatley it is not that easy. To draw a legend on the current plot, you need to run a new function called `legend()` and specify the appropriate arguments. The code to do so is provided below. Don't worry if it seems confusing, we plan on showing you a much more intuitive way of plotting your data.

```
legend("topleft", pch="*", col=c("blue", "green"), c("A", "B"), cex=0.8,
 	title="Celltype")
```
![scatter-5](../img/scatter-plot5.png)

***

**Exercise** 


1. Change the color scheme in the scatterplot, such that it reflects the `genotype` of samples rather than `celltype`.

2. Use R help to find out how to increase the size of the text on the axis labels.

***



## Advanced figures (`ggplot2`)


More recently, R users have moved away from base graphic options and towards a plotting package called [`ggplot2`](http://docs.ggplot2.org/). This package adds a lot of functionality to the basic plots described above. The syntax takes some getting used to, but once you have it you will find it's extremely powerful and flexible. We can start by re-creating the scatterplot but using ggplot functions to get a feel for the syntax.

`ggplot2` is best used on data in the `data.frame` form, so we will will work with `metadata` for the following figures. Let's start by loading the `ggplot2` library.

```{r}
library(ggplot2)
```

The `ggplot()` function is used to **initialize the basic graph structure**, then we add to it. The basic idea is that you specify different parts of the plot, and add them together using the `+` operator.

Let's start: 

```{r, eval=FALSE}
ggplot(new_metadata) # what happens? 
```
You get an blank plot, because you need to **specify layers** using the `+` operator.

One type of layer is **geometric objects**. These are the actual marks we put on a plot. Examples include:

* points (`geom_point`, for scatter plots, dot plots, etc)
* lines (`geom_line`, for time series, trend lines, etc)
* boxplot (`geom_boxplot`, for, well, boxplots!)

For a more exhaustive list on all possible geometric objects and when to use them check out [Hadley Wickham's RPubs](http://rpubs.com/hadley/ggplot2-layers). 

A plot **must have at least one geom**; there is no upper limit. You can add a geom to a plot using the `+` operator

```{r, eval=FALSE}
ggplot(new_metadata) +
  geom_point() # note what happens here
```

You will find that even though we have added a layer by specifying `geom_point`, we get an error. This is because each type of geom usually has a **required set of aesthetics** to be set. Aesthetic mappings are set with the aes() function and can be set inside `geom_point()` to be specifically applied to that layer. If we supplied aesthetics within `ggplot()`, they will be used as defaults for every layer. Examples of aesthetics include:

* position (i.e., on the x and y axes)
* color ("outside" color)
* fill ("inside" color) 
* shape (of points)
* linetype
* size

To start, we will add position for the x- and y-axis since `geom_point` requires the most basic information about a scatterplot, i.e. what you want to plot on the x and y axes. All of the others mentioned above are optional.

```{r, fig.align='center'}
ggplot(new_metadata) +
     geom_point(aes(x = age_in_days, y= samplemeans))
```

 ![ggscatter1](../img/ggscatter-1.png) 


Now that we have the required aesthetics, let's add some extras like color to the plot. We can `color` the points on the plot based on genotype, by specifying the column header. You will notice that there are a default set of colors that will be used so we do not have to specify. Also, the legend has been conveniently plotted for us!


```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype)) 
```

 ![ggscatter1.1](../img/ggscatter-2.png) 


Alternatively, we could color based on celltype by changing it to `color =celltype`. Let's try something different and have both **celltype and genotype identified on the plot**. To do this we can assign the `shape` aesthetic the column header, so that each celltype is plotted with a different shaped data point. Add in `shape = celltype` to your aesthetic and see how it changes your plot:

```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
  			shape=celltype)) 
`````

 ![ggscatter3](../img/ggscatter-3.png) 


The **size of the data points** are quite small. We can adjust that within the `geom_point()` layer, but does **not** need to be **included in `aes()`** since we are specifying how large we want the data points, rather than mapping it to a variable. Add in the `size` argument by specifying the magnitude relative to the default (`rel(3.0)`):

```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
  			shape=celltype), size=rel(3.0)) 
`````

 ![ggscatter4](../img/ggscatter-4.png)
  

The labels on the x- and y-axis are also quite small and hard to read. To change their size, we need to add an additional **theme layer**. The ggplot2 `theme` system handles non-data plot elements such as:

* Axis label aesthetics
* Plot background
* Facet label backround
* Legend appearance

There are built-in themes we can use (i.e. `theme_bw()`) that mostly change the background/foreground colours, by adding it as additional layer. Or we can adjust specific elements of the current default theme by adding the `theme()` layer and passing in arguments for the things we wish to change. Or we can use both.

Let's add a layer `theme_bw()`. Do the axis labels or the tick labels get any larger by changing themes? Not in this case. But we can add arguments using `theme()` to change it ourselves. Since we are adding this layer on top (i.e later in sequence), any features we change will override what is set in the `theme_bw()`. Here we'll increase the size of the axes labels and axes tick labels to be 1.5 times the default size.

```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
  			shape=celltype), size=rel(3.0)) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
  		axis.title = element_text(size=rel(1.5)))			
`````
 
 ![ggscatter5](../img/ggscatter-5.png)
 

***

**Exercise**

1. The current axis label text defaults to what we gave as input to `geom_point` (i.e the column headers). We can change this by adding additional layers called `xlab()` and `ylab()` for the x- and y-axis, respectively. Add these layers to the current plot such that the x-axis is labeled "Age (days)" and the y-axis is labeled "Mean expression".
2. Use the `ggtitle` layer to add a title to your plot.

***



## Histogram

To plot a histogram we require another type of geometric object called `geom_histogram`, which requires a statistical transformation. Some plot types (such as scatterplots) do not require transformations, each point is plotted at x and y coordinates equal to the original value. Other plots, such as boxplots, histograms, prediction lines etc. need to be transformed. Usually these objects have has a default statistic for the transformation, but that can be changed via the `stat_bin` argument. 

Let's plot a histogram of sample mean expression in our data:

```
ggplot(new_metadata) +
  geom_histogram(aes(x = samplemeans))
```


 ![gghist1](../img/gghist-1.png) 
 

You will notice that even though the histogram is plotted, R gives a warning message ``stat_bin()` using `bins = 30`. Pick better value with `binwidth`.` These are the transformations we discussed. Apparently the default is not good enough. 

Let's change the binwidth values. How does the plot differ?

```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_histogram(aes(x = samplemeans), stat = "bin", binwidth=0.8)
```

 ![gghist2](../img/gghist-2.png) 

***

**Exercise**

1. Use what you learned previously about the `theme()` layer to make the text larger for x- and y-axis (not the tick labels).
2. Also, add an appropriate title to this plot.

***



## Boxplot

Now that we have all the required information for plotting with ggplot2 let's try plotting a boxplot.

1. Use the `geom_boxplot()` layer to plot the differences in sample means between the Wt and KO genotypes.
2. Add a title to your plot.
3. Add 'Genotype' as your x-axis label and 'Mean expression' as your y-axis labels.
4. Change the size of your axes text to 1.5x larger than the default.
5. Change the size of your axes tick labels to 1.25x larger than the default.
6. Change the size of your plot title in the same way that you change the size of the axes text but use `plot.title`.

*BONUS: Use the `fill` aesthetic to look at differences in sample means between celltypes within each genotype.*


**Our final figure should look something like that provided below.**

 ![ggbox](../img/ggboxplot.png)


## Exporting figures to file

There are two ways in which figures and plots can be output to a file (rather than simply displaying on screen). The first (and easiest) is to export directly from the RStudio 'Plots' panel, by clicking on `Export` when the image is plotted. This will give you the option of `png` or `pdf` and selecting the directory to which you wish to save it to. 


The second option is to use R functions in the console, allowing you the flexibility to specify parameters to dictate the size and resolution of the output image, and allowing you to save multiple plots at once.  In R’s terminology, output is directed to a particular output device and that dictates the output format that will be produced.  A device must be created or “opened” in order to receive graphical output and, for devices that create a file
on disk, the device must also be closed in order to complete the output.

Let's print our scatterplot to a pdf file format. First you need to initialize a plot using a function which specifies the graphical format you intend on creating i.e.`pdf()`, `png()`, `tiff()` etc. Within the function you will need to specify a name for your image, and the with and height (optional). This will open up the device that you wish to write to:

	pdf("figure/scatterplot.pdf")


Then we plot the image to the device, using the ggplot scatterplot that we just created. 

```
ggplot(metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
  			shape=celltype), size=rel(3.0)) 
```

Finally, close the "device", or file, using the `dev.off()` function. There are also `bmp`, `tiff`, and `jpeg` functions, though the jpeg function has proven less stable than the others. 
  			
```    
dev.off()
```

***Note 1:*** *You will not be able to open and look at your file using standard methods (Adobe Acrobat or Preview etc.) until you execute the `dev.off()` function.*

***Note 2:*** *If you had made any additional plots before closing the device, they will all be stored in the same file; each plot usually gets its own page, unless you specify otherwise.*



---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
