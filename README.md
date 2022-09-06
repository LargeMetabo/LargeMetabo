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

    install_github("LargeMetabo/LargeMetabo", force = TRUE, build_vignettes = TRUE)
    library(LargeMetabo)

### Data Integration for Multiple Analytical Experiments

For data integration, multiple datasets from different analytical experiments can be used as the input of the LargeMetabo package. Before data integration, the csv files containing a feature-by-sample matrix should be prepared in advance. Each dataset (csv file) contains five essential columns providing the information of mass, retention time, intensity, isotope and adduct. The first two columns provide the mass and retention time, and samples must be kept in columns with the sample names in the first row.

    AlignData <- Integrate_Data(MutileGroup, RTTolerance1 = 10, mzTolerance1 = 0.1,
                               RTTolerance2 = 10, mzTolerance2 = 0.1)
    AlignData[1:5,1:5]

### Batch Effects Removal after Data Integration 

Various methods are provided in the LargeMetabo package for removing batch effects in different analytical experiments, including batch mean-centering (BMC/PAMR), the empirical Bayes method (ComBat/EB), and global normalization (GlobalNorm). An example input file for batch effect removal is provided in the LargeMetabo package.

    DataAfterBatch <- Removal_Batch(MutileAlign, n = 3, algorithm = "BMC/PAMR")
    DataAfterBatch[1:5,1:5]

### Sample Separation 

There are four sample separation methods for visualizing the clustering and separation of different samples. In the LargeMetabo package, the four methods are provided for sample separation. An example input file for sample separation is provided in the LargeMetabo package.

    finalData <- MarkerData$finalData
    finalLabel <- MarkerData$finalLabel
    Sample_Separation(finalData, finalLabel, clusters = 2, method = "HCA")

### Marker Identification 

In the marker identification step, there are 13 popular strategies to identify metabolic markers for the given dataset. These strategies include fold change (FC), partial least squares discrimination analysis (PLS-DA), orthogonal PLS-DA (OPLS-DA), Student’s t-test, Chi-squared test, correlation-based feature selection (CFS), entropy-based filter method, linear models and empirical Bayes method, recursive elimination of features (Relief), random forest-recursive feature elimination (RF-RFE), significance analysis for microarrays (SAM), support vector machine-recursive feature elimination (SVM-RFE), and Wilcoxon rank sum (WRS). An example input file for marker identification is provided in the LargeMetabo package.

    finalData <- MarkerData$finalData
    finalLabel <- MarkerData$finalLabel
    MarkerResult <- Marker_Identify(finalData, finalLabel, method = "FC")
    MarkerResult$FC_table[1:5,]

### Metabolite Annotation for MS1

When performing metabolite annotation for primary mass spectrometry (MS1), a compound list containing the studied m/z features should be properly provided. An example input for metabolite annotation for primary mass spectrometry is provided in the LargeMetabo package.

    AnnotaMS <- AnnotaData$AnnotaMS
    MetaboAResult <- Metabo_Annotation(AnnotaMS, masstole = 0.05, toleUnit = 1, annotaDB = "metlin",
                                   ionMode  = "pos")
    MetaboAResult$`M+H-2H2O`[1:5,]

### Metabolite Annotation for MS/MS

When performing metabolite annotation for tandem mass spectrometry (MS/MS), the information containing parent ion mass and MS/MS peak list (the first column is m/z value and the second column is the intensity) should be properly provided in this study. An example input for metabolite annotation for tandem mass spectrometry is provided in the LargeMetabo package. These example data embedded in the LargeMetabo package include the parent ion mass and MS/MS peak list (m/z & intensity).

    ParentMass <- AnnotaData$ParentMass
    TandemData <- AnnotaData$TandomData
    AnnotaParamTandem <- Annota_Tandem(ParentMass, TandemData, massTandem = 0.1, toleUnitTandem = 1,
                                   massmzTandem = 0.5, toleUnitmzTandem = 1, ModeTandem = "Positive",
                                   ionEnergy = "low(10V)")
    annota_Data_Tandem(AnnotaParamTandem)[1:5,]
    Annota_Tandem_plot(AnnotaParamTandem, TandemData)

### Enrichment Analysis for KEGG Pathways 

When performing enrichment analysis for KEGG pathways, a compound list should be properly provided. An example input for enrichment analysis for the KEGG pathways is provided in the LargeMetabo package.

    sampleDatakegg <- EnrichData$sampleDatakegg
    EnrichParam <- KEGG_Enrich_PlotPanel(sampleDatakegg, enrichDB = "kegg", pvalcutoff = 0.05,
                                   IDtype = 1, cateIdx = 1)
    EnrichResultList <- Enrichment(EnrichParam)
    EnrichFC <- seq(from = -2,to = 2, length.out = 24)
    KEGG_Enrich_Plot(EnrichResultList = EnrichResultList, cpdID = sampleDatakegg, cpdFC = EnrichFC)

### Enrichment Analysis for other Databases

When performing enrichment analysis for databases other than KEGG pathways, a compound list should also be properly provided. An example input for enrichment analysis for the classes of food components and food additives is provided in the LargeMetabo package.

    sampleDatacas <- EnrichData$sampleDatacas
    enrichDB <- EnrichData$enrichDB
    EnrichParam <- KEGG_Enrich_PlotPanel(sampleDatacas, enrichDB = enrichDB, pvalcutoff = 0.05,
                                   IDtype = 2, cateIdx = 1)
    EnrichResultList <- Enrichment(EnrichParam)
    dbChoice <- enrichDB
    Enrich_Plot(dbChoice, EnrichResultList)


