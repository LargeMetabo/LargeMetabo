# LargeMetabo: an out-of-the-box tool for processing and analyzing large-scale metabolomic data


### System requirements

Dependent on R (>= 3.5.0)

If you did not install the R software yet,you can download R >= 3.5.0  from https://www.r-project.org

### Installation
     
LargeMetabo package depends on several packages, which can be installed using the following commands in an R session:

    install.packages(c("corrplot", "e1071", "factoextra", "FSelector", "genefilter", "ggfortify", "ggplot2",
                       "igraph", "MASS", "mixOmics", "SOMbrero", "varSelRF"))
    
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

### Data Integration for Multiple Analytical Experiments

    AlignData <- Integrate_Data(MutileGroup, RTTolerance1 = 10, mzTolerance1 = 0.1,
                               RTTolerance2 = 10, mzTolerance2 = 0.1)
    AlignData[1:5,1:5]

### Batch Effects Removal after Data Integration 

    DataAfterBatch <- Removal_Batch(MutileAlign, n = 3, algorithm = "BMC/PAMR")
    DataAfterBatch[1:5,1:5]

### Sample Separation 

    finalData <- MarkerData$finalData
    finalLabel <- MarkerData$finalLabel
    Sample_Separation(finalData, finalLabel, clusters = 2, method = "HCA")

### Marker Identification 

    finalData <- MarkerData$finalData
    finalLabel <- MarkerData$finalLabel
    MarkerResult <- Marker_Identify(finalData, finalLabel, method = "FC")
    MarkerResult$FC_table[1:5,]

### Metabolite Annotation for MS1

    AnnotaMS <- AnnotaData$AnnotaMS
    MetaboAResult <- Metabo_Annotation(AnnotaMS, masstole = 10, toleUnit = 1, annotaDB = "metlin",
                                   ionMode  = "pos")
    MetaboAResult$`M+H-2H2O`[1:5,]

### Metabolite Annotation for MS/MS

    ParentMass <- AnnotaData$ParentMass
    TandemData <- AnnotaData$TandomData
    AnnotaParamTandem <- Annota_Tandem(ParentMass, TandemData, massTandem = 0.1, toleUnitTandem = 1,
                                   massmzTandem = 0.5, toleUnitmzTandem = 1, ModeTandem = "Positive",
                                   ionEnergy = "low(10V)")
    annota_Data_Tandem(AnnotaParamTandem)[1:5,]
    Annota_Tandem_plot(AnnotaParamTandem, TandemData)

### Enrichment Analysis for KEGG Pathways 

    sampleDatakegg <- EnrichData$sampleDatakegg
    EnrichParam <- KEGG_Enrich_PlotPanel(sampleDatakegg, enrichDB = "kegg", pvalcutoff = 0.05,
                                   IDtype = 1, cateIdx = 1)
    EnrichResultList <- Enrichment(EnrichParam)
    EnrichFC <- seq(from = -2,to = 2, length.out = 24)
    KEGG_Enrich_Plot(EnrichResultList = EnrichResultList, cpdID = sampleDatakegg, cpdFC = EnrichFC)

### Enrichment Analysis for other Databases

    sampleDatacas <- EnrichData$sampleDatacas
    enrichDB <- EnrichData$enrichDB
    EnrichParam <- KEGG_Enrich_PlotPanel(sampleDatacas, enrichDB = enrichDB, pvalcutoff = 0.05,
                                   IDtype = 2, cateIdx = 1)
    EnrichResultList <- Enrichment(EnrichParam)
    dbChoice <- enrichDB
    Enrich_Plot(dbChoice, EnrichResultList)


