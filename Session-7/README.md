# R Markdown and Quarto with RNA-seq

This week we're going to explore RMarkdown and its newer decendant Quarto. 

To get started, a nice introduction to RMarkdown can be found [here](https://rmarkdown.rstudio.com/articles_intro.html). 

A complete reference and guide can be gound [here](https://bookdown.org/yihui/rmarkdown/output-formats.html).

Additionally, [this PDF](https://www.rstudio.com/wp-content/uploads/2015/02/rmarkdown-cheatsheet.pdf) is a useful quick reference.

### Scenario

You are a scientist working with a team studying cancer-associated fibroblasts. 
One of your collaborators, a computational researcher, has just completed a bulk RNA-seq analysis, and has sent you the script they wrote for the initial part of the analysis. 
You want to show this analysis to two people, a data scientist on your team and the non-computational principle investigator of the research. 

We have two files, `analysis.R` and `analysis.docx`, which provide the analysis code and a small write-up of the analysis. 
However, we want to instead make a single integrated report we can send out. 

In order to show the analysis, you decide to create a report based on the script. 
However, you think that what is included in the report should probably be a bit different for your two audiences. 

Today, we're going to go through the steps of:

1. Converting R script to an RMarkdown file. 
2. Creating a RMarkdown report for our computational collaborator. 
3. Creating a RMarkdown report for our non-computational collaborator. 

### Experimental Summary

The paper this scenario is based on can be found here: [Induction of fibroblast senescence generates a non-fibrogenic myofibroblast phenotype that differentially impacts on cancer prognosis.](http://europepmc.org/article/MED/27992856)

The goal of this experiment was to compared the transcriptomics of myofibroblasts and senescent fibroblasts. To achieve this goal, RNA-seq experiment was carried out on human fetal foreskin fibroblasts 2 (HFFF2) treated with 2ng/ml TGF-beta-1 to induce myofibroblast differentiation or 10Gy gamma irradiation to induce senescence. RNA was isolated 7 days upon this treatments.

This report goes through some of the RNA-Seq analysis in R experiment starting with a read count matrix (processed using `salmon`). This report includes:

- Exploring count data after importing them into R and 
- Normalizing RNA-seq counts
- Conducting quality assessment of counts

## Part 1: Creating a report for a computational audience

0. Open and run analysis.R

  Let's make sure first that we have all the needed libraries install, we're in the correct working directory, and the analysis works. 
  After pulling the updated `Fall-2022` repository or downloading and extracting the repository's `.zip` file, open `analysis.R` and run all code. 

1. Create a new RMarkdown file
 
    Within R Studio, go to `File -> New File -> R Markdown`. 
    Leave `document` selected in the left sidebar of the pop-up window and `html` as the output format.
    Name the file `computational_report` and click `ok`.
  
2. Update header and setup block

    Your file should start with something like this right now: 

    ````
    ---
    title: "computational_report"
    output: html_document
    date: "2022-12-15"
    ---

    ```{r setup, include=FALSE}
    knitr::opts_chunk$set(echo = TRUE)
    ```
    ````

    Here, we have the header in [YAML format](https://zsmith27.github.io/rmarkdown_crash-course/lesson-4-yaml-headers.html) followed by the `setup` chunk. 

    We can now consider our global options. 
    While we want to output a `html_document` it can be nice to add other settings by including `html_notebook` in the output setting of the header. 
    We can then include options like adding a table of contents with `toc:true` and adding a theme. 

    We like the simple `cerulean` theme, but you can look at a list of other themes [here](https://www.datadreaming.org/post/r-markdown-theme-gallery/). 

    Looking at the reference and guide, what are some settings we might want globally for a report aimed at a computational audience?

    Settings to think about inside `knitr::opts_chunk$set` in the `setup` chunk are:

    - include = FALSE prevents code and results from appearing in the finished file. R Markdown still runs the code in the chunk, and the results can be used by other chunks.
    - echo = FALSE prevents code, but not the results from appearing in the finished file. This is a useful way to embed figures.
    - message = FALSE prevents messages that are generated by code from appearing in the finished file.
    - warning = FALSE prevents warnings that are generated by code from appearing in the finished.
    fig.cap = "..." adds a caption to graphical results.

    Do we want a computational collaborator to be able to see our code? 
    What about warnings? 

3. Bringing in the code

   For now, let's just paste the entire `analysis.R` script into our report as a single code block:

   ````
   ```{r full_analysis}
   library(DESeq2)
   # ... rest of code here
   ggbiplot::ggbiplot(vsd.prcp, scale = 0, var.axes = FALSE)
   ```
   ````
4. Adding date and session info
  
    Given that we are sending this report to a computational audience, we want to make sure we include detailed information about our packages and when we ran this analysis, in case our data scientist collaborator wants to re-run the analysis. 

    To do this, add the following text and separate chunk at the bottom of the markdown report:

    ````
    ## Wrap up

    Date stamp and session information. 

    ```{r date stamp}
    # add date stamp
    library(lubridate)
    lubridate::stamp("Data updated December 31, 1979")(lubridate::now())
    ```
    ```{r session info}
    # add session info
    library(utils)
    utils:::print.sessionInfo(sessionInfo()[-8]) 
    ```
    ````
  
5. [Draw the rest of the owl](https://knowyourmeme.com/memes/how-to-draw-an-owl)
  
   Now that we have our report setup, go though it and `analysis.docx`.
   Luckily your collaborator included some comments in their code, but we'll want to incorporate text from the analysis write-up to provide more context. 

   Consider how you to want to split up your code into meaningful chunks. 
   Are there particular parts of the code where you would want to change the global settings you chose?
   
   It would be nice in a few places to show some tables of the data. 
   You can accomplish this using the `knitr::kable` function. 
   For instance, adding the lines:

   ```
   knitr::kable(head(counts, 10), caption = "First 10 rows of raw count data")
   knitr::kable(sampleinfo, caption = "Sample Metadata")
   ```
   Would allow you to preview your counts and sample information tables after you load them. 
   
   For blocks which plot figures, block-level settings such as `fig.cap` can be used to adjust and add to plots. 
   
   We also currently try PCA in two different ways, using `PCAPlot` and using `prcomp` with only the 500 top genes. 
   One nice way to show these kinds of multiple but parallel options is to [use tabs in the report](https://bookdown.org/yihui/rmarkdown-cookbook/html-tabs.html). 


   __Remember to periodicly try knitting your code report!__


## Part 2: Creating a report for a non-computational audience

We are now going to create another report, this time for a non-computational audience. 

This collaborator, the PI of the project, isn't as interested in experimental details or reproducing the experiment themselves. 
They simply want to know what you did and what the results are. 

Create a new R Markdown file called `biological_report`, copying your comptuational report. 

Consider the settings at the top of the file, is there anything we want to change?
Go through your existing text and output, what might we want to add or remove for this report?

Try adding at least 1 line where you use [inline code](https://bookdown.org/yihui/rmarkdown-cookbook/r-code.html) to create a dynamic sentence based on the data, such as generating a sentence that says:

> Our dataset includes counts of `[number of genes]` genes for `[number of samples]` samples. 
