# LargeMetabo: an out-of-the-box tool for processing and analyzing large-scale metabolomic data


### System requirements

Dependent on R (>= 3.5.0)

If you did not install the R software yet,you can download R >= 3.5.0  from https://www.r-project.org

### Installation
     
LargeMetabo package depends on several packages, which can be installed using the following commands in an R session:

    install.packages(c("corrplot", "e1071", "factoextra", "FSelector", "genefilter", "ggfortify", "ggplot2", "igraph", 
                       "MASS", "mixOmics", "SOMbrero", "varSelRF"))
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c("CluMSID","genefilter","ropls","siggenes", "GenomeInfoDbData"))
    
    if (!require("devtools")) install.packages("devtools")
    devtools::install_github("rstudio/d3heatmap")

The LargeMetabo package is provided through GitHub. In order to install it, devtools package available in CRAN (https://cran.r-project.org/) is required. To install devtools, the user must type the following commands in an R session:
    
    install.packages("devtools")
    library(devtools)

Once devtools package has been installed, the user can install LargeMetabo package by typing the following commands in an R session:

    install_github("LargeMetabo/LargeMetabo", force = TRUE)
    library(LargeMetabo)
