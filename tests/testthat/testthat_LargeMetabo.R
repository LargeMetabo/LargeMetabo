library(LargeMetabo)
context("LargeMetabo: an out-of-the-box tool for processing and analyzing large-scale metabolomic data")





#### Biomarker identification

markers <- Marker_Assess(t(iris[1:100,1:4]), as.numeric(iris[1:100,5]), method = "PLS-DA")

test_that("Biomarker identification", {
  expect_equal(markers, "The AUC of the ROC curve is 1.000")
})




#### Metabolite annotation

ParentMass <- 181
TandemData <- matrix(c(122,123,140,141,142,182,183,184,0.59,2.97,100,2.71,2.72,36.24,2.38,1.17), 8, 2)
colnames(TandemData) <- c("V1", "V2")
row.names(TandemData) <- 1:8


AnnotaParamTandem <- Annota_Tandem(ParentMass, TandemData)


test_that("Metabolite annotation", {
  expect_equal(annota_Data_Tandem(AnnotaParamTandem)[1,5], 1)
})



#### Enrichment analysis

sampleDatakegg <- c("C00001","C00002","C00003","C00004","C00005","C00006","C00007","C00008","C00009","C00010","C00011","C00012","C00013","C00014","C00015","C00016","C00017",
                    "C00018","C00019","C00020","C00021","C00022","C00023","C00024","C00025","C00026","C00027","C00028","C00029","C00030","C00031","C00032","C00033","C00034",
                    "C00035","C00036","C00037","C00038","C00039","C00040","C00041","C00042","C00043","C00044","C00045","C00046","C00047","C00048","C00049","C00050","C00255",
                    "C00256","C00257","C00258","C00259","C00261","C00262","C00263","C00264","C00265","C00266","C00267","C00268","C00269","C00270","C00272","C00273","C00275",
                    "C00279","C00280","C00282","C00283","C00284","C00286","C00288","C00290","C00291","C00292","C00294","C00295","C00296","C00297","C00298","C00299","C00300",
                    "C00301","C00302","C00303","C00304","C00305","C00306","C00307","C00308","C00309","C00310","C00311","C00312","C00313","C00314","C00315","C00317")


EnrichParam <- KEGG_Enrich_PlotPanel(sampleDatakegg, enrichDB = "kegg", pvalcutoff = 0.05, IDtype = 1, cateIdx = 1)

EnrichResultList <- Enrichment(EnrichParam)


test_that("Enrichment analysis", {
  expect_equal(EnrichResultList$Table.Result[1,1], "hsa00630")
})













