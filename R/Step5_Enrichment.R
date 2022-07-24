
#' Title KEGGEnrichPlotPanel
#'
#' @param sampleData The charactor of input for metabolite enrichment.
#' @param enrichDB The database of the input data. The database could be kegg, smpdb, cfam, foodb, biofunc, tcm, spectax or toxin.
#' @param pvalcutoff The cutoff of p value for metabolite enrichment.
#' @param IDtype The number for name type for metabolite enrichment, such as KEGG ID, CAS ID, PubChem ID, Name or HMDB ID.
#' @param cateIdx The number of category for metabolite enrichment.
#'
#' @return The parameters for metabolite enrichment.
#' @export
#'
#' @examples sampleDatakegg <- EnrichData$sampleDatakegg
#' KEGGEnrichPlotPanel(sampleDatakegg, enrichDB = "kegg", pvalcutoff = 0.05, IDtype = 1, cateIdx = 1)

KEGGEnrichPlotPanel <- function(sampleData, enrichDB = "kegg", pvalcutoff = 0.05, IDtype = 1, cateIdx = 1){


  #KEGGEnrichPlotPanel

  library(ggplot2)


  param <- list()

  ### padjIdx == 1,2,3,4,5,6,7,8
  padjIdx <- 1
  PAdjMethodList <- c("bonferroni","hochberg","hommel","holm","BH","BY","fdr","none")
  param$padj <- PAdjMethodList[padjIdx]
  param$pvaladjmethod <- PAdjMethodList[padjIdx]

  ## enrichDB==kegg, smpdb, cfam, foodb, biofunc, tcm, spectax, toxin
  enrichDB <- enrichDB
  # enrichDB <- "toxin"
  param$db <- as.character(enrichDB)

  ## pvalcutoff==0.05
  pvalcutoff <- pvalcutoff
  param$pvalcutoff <- as.numeric(pvalcutoff)

  ## sampleData==输入的化合物名称
  sampleData <- sampleData
  # sampleData <- read.table("clipboard",header=FALSE)[,1]
  EnrichCompound <- unique(sampleData)

  ## inputtype
  if (enrichDB=='kegg'){
    tochoice <- c("KEGG_ID","CAS_ID","PubChem_ID","Name")
  }else if (enrichDB=='smpdb'){
    tochoice <- c("HMDB_ID", "KEGG_ID","CAS_ID","PubChem_ID","Name", "Pw_ID")
  }else if (enrichDB=='cfam'){
    tochoice <- c("CFAM_ID","CAS_ID","PubChem_ID","Name", "HMDB_ID")
  }else if (enrichDB=='foodb'){
    tochoice <- c("Food_Agr_ID","CAS_ID","PubChem_ID","Name", "HMDB_ID")
  }else if (enrichDB=='biofunc'){
    tochoice <- c("HMDB_ID","KEGG_ID","CAS_ID","PubChem_ID","Name")
  }else if (enrichDB=='tcm'){
    tochoice <- c("MMEASE_ID","PubChem_ID","Name", "CAS_ID")
  }else if (enrichDB=='spectax'){
    tochoice <- c("MMEASE_ID","PubChem_ID","Name", "CAS_ID","HMDB_ID")
  }else if (enrichDB=='toxin'){
    tochoice <- c("T3DB_ID","PubChem_ID","Name", "CAS_ID","HMDB_ID","KEGG_ID")
  }

  #inputtype <- tochoice[1]
  inputtype <- tochoice[as.numeric(IDtype)]

  ## IDTrans

  db <- param$db
  idtype <- inputtype
  cpdlist <- EnrichCompound



  IDTrans <- function(db,idtype,cpdlist){

    #load("D:/Code/Rpackage/LargeMetabo/data/tableName_all.RData")
    
    
    
    tableName_kegg <- function() {
      utils::data(list="tableName_kegg", package="LargeMetabo")
      get("tableName_kegg", envir = .GlobalEnv)
    }
    tableName_kegg <- tableName_kegg()
    
    
    tableName_smpdb <- function() {
      utils::data(list="tableName_smpdb", package="LargeMetabo")
      get("tableName_smpdb", envir = .GlobalEnv)
    }
    tableName_smpdb <- tableName_smpdb()
    
    
    tableName_cfam <- function() {
      utils::data(list="tableName_cfam", package="LargeMetabo")
      get("tableName_cfam", envir = .GlobalEnv)
    }
    tableName_cfam <- tableName_cfam()
    
    
    
    tableName_biofunc <- function() {
      utils::data(list="tableName_biofunc", package="LargeMetabo")
      get("tableName_biofunc", envir = .GlobalEnv)
    }
    tableName_biofunc <- tableName_biofunc()
    
    
    tableName_foodb <- function() {
      utils::data(list="tableName_foodb", package="LargeMetabo")
      get("tableName_foodb", envir = .GlobalEnv)
    }
    tableName_foodb <- tableName_foodb()
    
    
    tableName_tcm <- function() {
      utils::data(list="tableName_tcm", package="LargeMetabo")
      get("tableName_tcm", envir = .GlobalEnv)
    }
    tableName_tcm <- tableName_tcm()
    
    
    tableName_spectax <- function() {
      utils::data(list="tableName_spectax", package="LargeMetabo")
      get("tableName_spectax", envir = .GlobalEnv)
    }
    tableName_spectax <- tableName_spectax()
    
    
    tableName_toxin <- function() {
      utils::data(list="tableName_toxin", package="LargeMetabo")
      get("tableName_toxin", envir = .GlobalEnv)
    }
    tableName_toxin <- tableName_toxin()
    
    
    
    tableName_all <- list()
    
    tableName_all$tableName_kegg <- tableName_kegg
    tableName_all$tableName_smpdb <- tableName_smpdb
    tableName_all$tableName_cfam <- tableName_cfam
    tableName_all$tableName_biofunc <- tableName_biofunc
    tableName_all$tableName_foodb <- tableName_foodb
    tableName_all$tableName_tcm <- tableName_tcm
    tableName_all$tableName_spectax <- tableName_spectax
    tableName_all$tableName_toxin <- tableName_toxin
    
    


    destColName <- NULL
    tableName <- NULL
    if (db=="kegg"){
      destColName <- "KEGG"
      tableName <- tableName_all$tableName_kegg

      cpds <- c()

      for (sourceID in cpdlist){
        if (idtype=="Name"){
          IDsql <- tableName$KEGG_ID[grep(sourceID, tableName$Name)]
        }else{
          IDsql <- tableName$KEGG_ID[grep(sourceID, tableName[,idtype])]
        }
        if (length(IDsql)>0){
          cpds <- c(cpds, as.character(IDsql))
        }
      }

    }else if (db=="smpdb"){
      destColName <- "HMDB"
      tableName <- tableName_all$tableName_smpdb

      cpds <- c()

      for (sourceID in cpdlist){
        if (idtype=="Name"){
          IDsql <- tableName$HMDB_ID[grep(sourceID, tableName$Name)]
        }else{
          IDsql <- tableName$HMDB_ID[grep(sourceID, tableName[,idtype])]
        }
        if (length(IDsql)>0){
          cpds <- c(cpds, as.character(IDsql))
        }
      }

    }else if (db=="cfam"){
      destColName <- "CFAM"
      tableName <- tableName_all$tableName_cfam

      cpds <- c()

      for (sourceID in cpdlist){
        if (idtype=="Name"){
          IDsql <- tableName$CFAM_ID[grep(sourceID, tableName$Name)]
        }else{
          IDsql <- tableName$CFAM_ID[grep(sourceID, tableName[,idtype])]
        }
        if (length(IDsql)>0){
          cpds <- c(cpds, as.character(IDsql))
        }
      }

    }else if (db=="biofunc"){
      destColName <- "HMDB"
      tableName <- tableName_all$tableName_biofunc

      cpds <- c()

      for (sourceID in cpdlist){
        if (idtype=="Name"){
          IDsql <- tableName$HMDB_ID[grep(sourceID, tableName$Name)]
        }else{
          IDsql <- tableName$HMDB_ID[grep(sourceID, tableName[,idtype])]
        }
        if (length(IDsql)>0){
          cpds <- c(cpds, as.character(IDsql))
        }
      }

    }else if (db=="foodb"){
      destColName <- "FoodAgr"
      tableName <- tableName_all$tableName_foodb

      cpds <- c()

      for (sourceID in cpdlist){
        if (idtype=="Name"){
          IDsql <- tableName$Food_Agr_ID[grep(sourceID, tableName$Name)]
        }else{
          IDsql <- tableName$Food_Agr_ID[grep(sourceID, tableName[,idtype])]
        }
        if (length(IDsql)>0){
          cpds <- c(cpds, as.character(IDsql))
        }
      }

    }else if (db=="tcm"){
      destColName <- "MMEASE"
      tableName <- tableName_all$tableName_tcm

      cpds <- c()

      for (sourceID in cpdlist){
        if (idtype=="Name"){
          IDsql <- tableName$MMEASE_ID[grep(sourceID, tableName$Name)]
        }else{
          IDsql <- tableName$MMEASE_ID[grep(sourceID, tableName[,idtype])]
        }
        if (length(IDsql)>0){
          cpds <- c(cpds, as.character(IDsql))
        }
      }

    }else if (db=="spectax"){
      destColName <- "MMEASE"
      tableName <- tableName_all$tableName_spectax

      cpds <- c()

      for (sourceID in cpdlist){
        if (idtype=="Name"){
          IDsql <- tableName$MMEASE_ID[grep(sourceID, tableName$Name)]
        }else{
          IDsql <- tableName$MMEASE_ID[grep(sourceID, tableName[,idtype])]
        }
        if (length(IDsql)>0){
          cpds <- c(cpds, as.character(IDsql))
        }
      }

    }else if (db=="toxin"){
      destColName <- "T3DB"
      tableName <- tableName_all$tableName_toxin

      cpds <- c()

      for (sourceID in cpdlist){
        if (idtype=="Name"){
          IDsql <- tableName$T3DB_ID[grep(sourceID, tableName$Name)]
        }else{
          IDsql <- tableName$T3DB_ID[grep(sourceID, tableName[,idtype])]
        }
        if (length(IDsql)>0){
          cpds <- c(cpds, as.character(IDsql))
        }
      }

    }

    return (cpds)
  }


  ### cpdlist
  param$cpdlist <- IDTrans(param$db,inputtype,EnrichCompound)




  if (enrichDB=="kegg"){
    SpeciesList <- c("hsa", "mmu", "rno", "bta", "gga", "dre", "dme", "cel", "sce", "osa", "ath",
                     "smm", "pfa", "tbr", "eco", "bsu", "ppu", "sau", "tma", "syf", "mlo")
    ## KEGGSpec=="Homo sapiens (human)"=1, "Mus musculus (mouse)"=2,
    ## "Rattus norvegicus (rat)"=3,"Bos taurus (cow)"=4,
    ## "Gallus gallus (chicken)"=5,"Danio rerio (zebrafish)"=6,
    ## "Drosophila melanogaster (fruit fly)"=7,"Caenorhabditis elegans (nematode)"=8,
    ## "Saccharomyces cerevisiae (yeast)"=9,"Oryza sativa japonica (Japanese rice)"=10,
    ## "Arabidopsis thaliana (thale cress)"=11,"Schistosoma mansoni"=12,
    ## "Plasmodium falciparum 3D7 (Malaria)"=13,"Trypanosoma brucei" =14,
    ## "Escherichia coli K-12 MG1655"=15,"Bacillus subtilis"=16,"Pseudomonas putida KT2440" =17,
    ## "Staphylococcus aureus N315(MRSA/VSSA)"=18,"Thermotoga maritima"=19,
    ## "Synechococcus elongatus PCC7942"=20,"Mesorhizobium loti" =21
    KEGGSpec <- 1
    specIdx <- as.numeric(KEGGSpec)
    param$organism <- SpeciesList[specIdx]

  }else if (enrichDB=="cfam"){
    CategoryList <-c("Class", "Superfamily", "Family")
    ## CFamCategory=="Class"=1, "Superfamily"=2, "Family"=3
    CFamCategory <- 1
    cateIdx <- as.numeric(cateIdx)
    param$cfam_category <- CategoryList[cateIdx]

  }else if (enrichDB=="biofunc"){
    CategoryList <-c("Blood","Urine","Cerebrospinal")
    ## Biofluid=="Blood"=1,"Urine"=2,"Cerebrospinal"=3
    Biofluid <- 1
    cateIdx <- as.numeric(cateIdx)
    param$biofluid <- CategoryList[cateIdx]

  }else if (enrichDB=="foodb"){
    CategoryList <-c("Food group","Food subgroup")
    ## foodb_category_type=="Food group"=1,"Food subgroup"=2
    foodb_category_type <- 1
    cateIdx <- as.numeric(cateIdx)
    param$category_type <- CategoryList[cateIdx]

  }else if (enrichDB=="tcm"){
    CategoryList <-c("TCM primary category","TCM secondary category")
    ## TCMClass=="TCM primary category"=1,"TCM secondary category"=2
    TCMClass <- 1
    cateIdx <- as.numeric(cateIdx)
    param$category_type <- CategoryList[cateIdx]

  }else if (enrichDB=="spectax"){
    CategoryList <-c("Superkingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    ## spectax_category_type=="Superkingdom"=1,"Kingdom"=2,"Phylum"=3,"Class"=4,"Order"=5,"Family"=6,"Genus"=7,"Species"=8
    spectax_category_type <- 1
    cateIdx <- as.numeric(cateIdx)
    param$category_type <- CategoryList[cateIdx]
  }

  #####
  EnrichParam <- param

  return(EnrichParam)


}







######### table #################

#' Title
#'
#' @param param The parameters for metabolite enrichment.
#'
#' @return The table of results for metabolite enrichment.
#' @export
#'
#' @examples sampleDatakegg <- EnrichData$sampleDatakegg
#' EnrichParam <- KEGGEnrichPlotPanel(sampleDatakegg, enrichDB = "kegg", pvalcutoff = 0.05, IDtype = 1, cateIdx = 1)
#' Enrichment(EnrichParam)

Enrichment <- function(param){

  #load("D:/Code/Rpackage/LargeMetabo/data/sql_all.RData")


  sql_all <- function() {
    utils::data(list="sql_all", package="LargeMetabo")
    get("sql_all", envir = .GlobalEnv)
  }

  sql_all <- sql_all()



  if (param$db=="kegg"){
    sql <- sql_all$sql_kegg
    namesql <- sql$pathway_name[grep(param$organism, sql$Orgnism)]
    CategoryName <- namesql

    backsql <- sql$Compound_ID[grep(param$organism, sql$Orgnism)]
    CategoryCpds <- as.list(backsql)
    CategoryCpds <- lapply(CategoryCpds,function(x){ strsplit(x, ";")[[1]]})

    idsql <- sql$Pathway_ID[grep(param$organism, sql$Orgnism)]
    names(CategoryCpds) <- idsql

  }else if (param$db=="smpdb"){
    sql <- sql_all$sql_smpdb
    namesql <- sql$SMPDB_Name
    CategoryName <- namesql

    backsql <- sql$HMDB_ID
    CategoryCpds <- as.list(backsql)
    CategoryCpds <- lapply(CategoryCpds,function(x){ strsplit(x, ";")[[1]]})

    idsql <- sql$SMPDB_ID
    names(CategoryCpds) <- idsql

  }else if (param$db=="cfam"){
    sql <- sql_all$sql_cfam
    namesql <- sql$Category_name[grep(param$cfam_category, sql$CFam_Class)]
    CategoryName <- namesql

    backsql <- sql$Cfam_Mol_ID[grep(param$cfam_category, sql$CFam_Class)]
    CategoryCpds <- as.list(backsql)
    CategoryCpds <- lapply(CategoryCpds,function(x){ strsplit(x, ";")[[1]]})

    idsql <- sql$Category_ID[grep(param$cfam_category, sql$CFam_Class)]
    names(CategoryCpds) <- idsql

  }else if (param$db=="biofunc"){
    sql <- sql_all$sql_biofunc
    namesql <- sql$Biofunction_name[grep(param$biofluid, sql$Biofluid)]
    CategoryName <- namesql

    backsql <- sql$HMDB_ID[grep(param$biofluid, sql$Biofluid)]
    CategoryCpds <- as.list(backsql)
    CategoryCpds <- lapply(CategoryCpds,function(x){ strsplit(x, ";")[[1]]})

    idsql <- sql$Biofunction_name[grep(param$biofluid, sql$Biofluid)]
    names(CategoryCpds) <- idsql

  }else if(param$db=="foodb"){
    sql <- sql_all$sql_foodb
    namesql <- sql$category_name[grep(param$category_type, sql$category_type)]
    CategoryName <- namesql

    backsql <- sql$cpd_id[grep(param$category_type, sql$category_type)]
    CategoryCpds <- as.list(backsql)
    CategoryCpds <- lapply(CategoryCpds,function(x){ strsplit(x, ";")[[1]]})

    idsql <- sql$category_name[grep(param$category_type, sql$category_type)]
    names(CategoryCpds) <- idsql

  }else if(param$db=="tcm"){
    sql <- sql_all$sql_tcm
    namesql <- sql$category_name[grep(param$category_type, sql$category_type)]
    CategoryName <- namesql

    backsql <- sql$cpd_id[grep(param$category_type, sql$category_type)]
    CategoryCpds <- as.list(backsql)
    CategoryCpds <- lapply(CategoryCpds,function(x){ strsplit(x, ";")[[1]]})

    idsql <- sql$category_name[grep(param$category_type, sql$category_type)]
    names(CategoryCpds) <- idsql

  }else if(param$db=="spectax"){
    sql <- sql_all$sql_spectax
    namesql <- sql$category_name[grep(param$category_type, sql$category_type)]
    CategoryName <- namesql

    backsql <- sql$cpd_id[grep(param$category_type, sql$category_type)]
    CategoryCpds <- as.list(backsql)
    CategoryCpds <- lapply(CategoryCpds,function(x){ strsplit(x, ";")[[1]]})

    idsql <- sql$category_name[grep(param$category_type, sql$category_type)]
    names(CategoryCpds) <- idsql

  }else if (param$db=="toxin"){
    sql <- sql_all$sql_toxin
    namesql <- sql$category_name
    CategoryName <- namesql

    backsql <- sql$cpd_id
    CategoryCpds <- as.list(backsql)
    CategoryCpds <- lapply(CategoryCpds,function(x){ strsplit(x, ";")[[1]]})

    idsql <- sql$category_name
    names(CategoryCpds) <- idsql

  }


  BackgroundCpds <- (unique(unlist(CategoryCpds)))

  ## Matched and unmatched Cpds
  MatchCpds <- intersect(param$cpdlist, BackgroundCpds)
  UnmatchCpds <- setdiff(param$cpdlist, MatchCpds)

  ## QueryCpds
  QueryCpds <- lapply(CategoryCpds, function(x){intersect(x,MatchCpds)})
  NonNullLines <- unlist(lapply(QueryCpds, function(x){ifelse(length(x)==0,F,T)}  ))
  QueryCpds <- QueryCpds[NonNullLines]

  ## CategroyCpds
  QueryPos <- match( names(QueryCpds), names(CategoryCpds) )
  CategoryCpds <-CategoryCpds[QueryPos]

  ## HyperGeometric Parameter
  #basic parameter
  x <- sapply(QueryCpds, length)
  m <- sapply(CategoryCpds, length)
  n <- length(BackgroundCpds) - m
  k <- rep(length(MatchCpds), times = length(m))
  args.df <- data.frame(x = x,m = m,n = n,k = k)
  #P-Value
  getPvalue <- function(obj) {phyper(obj[1] - 1, obj[2], obj[3], obj[4], lower.tail = FALSE)}
  Pvalues   <- apply(args.df, 1, getPvalue)
  PvaluesAdj<- p.adjust(Pvalues, method = param$pvaladjmethod)
  #Other pamrameter
  CpdsRatio <- apply(data.frame(a = x, b = k), 1, function(x){paste(x[1], "/", x[2], sep="", collapse="")})
  BackgroundRatio  <- apply(data.frame(a = m, b = m + n), 1, function(x){paste(x[1], "/", x[2], sep="", collapse="")})
  CategoryID <- names(QueryCpds)

  if(param$db=="kegg"){
    CategoryID <- gsub(pattern="map", replacement=param$organism, CategoryID)
  }
  CategoryName <- CategoryName[QueryPos]
  CompoundID <- sapply(QueryCpds, function(i) paste(i, collapse=" "))
  Count <- sapply(x,  function(x){x[1]})

  ## Build a result dataframe
  #----------------------------------------------------------------------------------------+
  if (param$db=="kegg"){
    EnrichResult <- data.frame(Pathway.ID = CategoryID,Pathway.name = CategoryName,
                               Count = Count, Mol.Ratio = CpdsRatio,Bg.Ratio  = BackgroundRatio,
                               P.Val = Pvalues,P.Val.adj = PvaluesAdj, Compound.ID = CompoundID)
  }else if (param$db=="smpdb"){
    EnrichResult <- data.frame(SMPBD.ID = CategoryID, SMPBD.name = CategoryName,
                               Count = Count, Mol.Ratio = CpdsRatio,Bg.Ratio  = BackgroundRatio,
                               P.Val = Pvalues,P.Val.adj = PvaluesAdj, Compound.ID = CompoundID)
  }else if (param$db=="cfam"){
    EnrichResult <- data.frame(Category.ID = CategoryID,Category.name = CategoryName,
                               Count = Count, Mol.Ratio = CpdsRatio,Bg.Ratio  = BackgroundRatio,
                               P.Val = Pvalues,P.Val.adj = PvaluesAdj, CFam.Mol.ID = CompoundID)
  }else if (param$db=="biofunc"){
    EnrichResult <- data.frame(Biofunction.name = CategoryName,
                               Count = Count, Mol.Ratio = CpdsRatio,Bg.Ratio  = BackgroundRatio,
                               P.Val = Pvalues,P.Val.adj = PvaluesAdj, HMDB.ID = CompoundID)
  }else if (param$db=="foodb"){
    EnrichResult <- data.frame(Category.name = CategoryName,
                               Count = Count, Mol.Ratio = CpdsRatio,Bg.Ratio  = BackgroundRatio,
                               P.Val = Pvalues,P.Val.adj = PvaluesAdj, HMDB.ID = CompoundID)
  }else if (param$db=="tcm"){
    EnrichResult <- data.frame(Category.name = CategoryName,
                               Count = Count, Mol.Ratio = CpdsRatio,Bg.Ratio  = BackgroundRatio,
                               P.Val = Pvalues,P.Val.adj = PvaluesAdj, HMDB.ID = CompoundID)
  }else if (param$db=="spectax"){
    EnrichResult <- data.frame(Category.name = CategoryName,
                               Count = Count, Mol.Ratio = CpdsRatio,Bg.Ratio  = BackgroundRatio,
                               P.Val = Pvalues,P.Val.adj = PvaluesAdj, HMDB.ID = CompoundID)
  }else if (param$db=="toxin"){
    EnrichResult <- data.frame(Category.name = CategoryName,
                               Count = Count, Mol.Ratio = CpdsRatio,Bg.Ratio  = BackgroundRatio,
                               P.Val = Pvalues,P.Val.adj = PvaluesAdj, Compound.ID = CompoundID)
  }
  #----------------------------------------------------------------------------------------|
  rownames(EnrichResult) <- NULL
  EnrichResult <- EnrichResult[ order(EnrichResult$P.Val),  ]

  ## Filter and Sort on result dataframe
  EnrichResult       <- EnrichResult[EnrichResult$P.Val <= param$pvalcutoff, ]
  EnrichResult$P.Val <- format(EnrichResult$P.Val, digits = 3, scientific = TRUE)
  EnrichResult$P.Val.adj <- format(EnrichResult$P.Val.adj, digits = 3, scientific = TRUE)

  ## Returning a result list with three parts.
  EnrichResultList <- list(Table.Result = EnrichResult, Match.cpd = MatchCpds, Unmatch.Mol = UnmatchCpds)


  return (EnrichResultList)

}
# Enrichment(param)





######### plot #################

#' Title KEGGEnrichPlot
#'
#' @param EnrichResultList The table of results for metabolite enrichment using KEGG database.
#' @param cpdID The charactor of input for metabolite enrichment.
#' @param cpdFC The distance in the metabolite enrichment plot using KEGG database.
#'
#' @return The plot of results for metabolite enrichment using KEGG database.
#' @export
#'
#' @examples sampleDatakegg <- EnrichData$sampleDatakegg
#' EnrichParam <- KEGGEnrichPlotPanel(sampleDatakegg, enrichDB = "kegg", pvalcutoff = 0.05, IDtype = 1, cateIdx = 1)
#' EnrichResultList <- Enrichment(EnrichParam)
#' EnrichFC <- seq(from=-2,to=2, length.out =24)
#' KEGGEnrichPlot(EnrichResultList=EnrichResultList,cpdID=sampleDatakegg,cpdFC=EnrichFC)

KEGGEnrichPlot <- function(EnrichResultList, cpdID, cpdFC) {
  require(ggplot2)

  s =2
  e = 10

  rn <- length(EnrichResultList$Match.cpd)
  cn <- nrow(EnrichResultList$Table.Result)
  # Constructing a matrix, containing rn rows and cn columns.
  datforplot <- matrix(0, rn, cn)
  rownames(datforplot) <- EnrichResultList$Match.cpd
  colnames(datforplot) <- as.character(EnrichResultList$Table.Result$Pathway.name)
  match.f <- function(x, y) {
    ifelse(length(grep(x, y)) != 1, 0, 1)
  }

  # Matched compounds.
  comps <- EnrichResultList$Match.cpd
  path.term <- EnrichResultList$Table.Result$Compound.ID
  for (i in 1:nrow(datforplot)) {
    for (path in 1:ncol(datforplot)) {
      datforplot[i, path] <- match.f(comps[i], path.term[path])
    }
  }

  # datforplot
  # Selecting several terms to visualize, such as from 1 to 8.
  ect.f <- function(x) {if (mean(x) == 0) FALSE else TRUE}
  dat.vis <- function(dat, s = 1, t = 8) {
    tmp <- dat[, s:t]
    tmp[apply(tmp, 1, ect.f), ]
  }

  dat.chord <- dat.vis(datforplot, s, e)
  dat.chord <- as.data.frame(dat.chord)
  dat.chord$logFC <- cpdFC[match(rownames(dat.chord),  as.character(cpdID))]# seq(-2, 2, length = nrow(dat.chord))

  #source("C:/Users/Yang/Documents/mmease001/R/GOCore.R")

  theme_blank <- theme(axis.line = element_blank(), axis.text.x = element_blank(),
                       axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
                       axis.title.y = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())



  draw_table <- function(data, col = ''){
    id <- term <- NULL
    colnames(data) <- tolower(colnames(data))
    if (length(col) == 1){
      tt1 <- ttheme_default()
    }else{
      text.col <- c(rep(col[1], sum(data$category == 'BP')), rep(col[2], sum(data$category == 'CC')), rep(col[3], sum(data$category == 'MF')))
      tt1 <- ttheme_minimal(core=list(fg.par = list(size = 4), bg.par = list(fill = text.col, col=NA, alpha= 1/3)), colhead=list(fg.par=list(col="black")))
    }
    table <- tableGrob(subset(data, select = c(id, term)), cols = c('ID', 'Describtion'), rows = NULL, theme = tt1)
    return(table)
  }



  bezier <- function(data, process.col){
    x <- c()
    y <- c()
    Id <- c()
    sequ <- seq(0, 1, by = 0.01)
    N <- dim(data)[1]
    sN <- seq(1, N, by = 2)
    if (process.col[1] == '') col_rain <- grDevices::rainbow(N) else col_rain <- process.col
    for (n in sN){
      xval <- c(); xval2 <- c(); yval <- c(); yval2 <- c()
      for (t in sequ){
        xva <- (1 - t) * (1 - t) * data$x.start[n] + t * t * data$x.end[n]
        xval <- c(xval, xva)
        xva2 <- (1 - t) * (1 - t) * data$x.start[n + 1] + t * t * data$x.end[n + 1]
        xval2 <- c(xval2, xva2)
        yva <- (1 - t) * (1 - t) * data$y.start[n] + t * t * data$y.end[n]
        yval <- c(yval, yva)
        yva2 <- (1 - t) * (1 - t) * data$y.start[n + 1] + t * t * data$y.end[n + 1]
        yval2 <- c(yval2, yva2)
      }
      x <- c(x, xval, rev(xval2))
      y <- c(y, yval, rev(yval2))
      Id <- c(Id, rep(n, 2 * length(sequ)))
    }
    df <- data.frame(lx = x, ly = y, ID = Id)
    return(df)
  }



  circle_dat <- function(terms, genes){

    colnames(terms) <- tolower(colnames(terms))
    terms$genes <- toupper(terms$genes)
    genes$ID <- toupper(genes$ID)
    tgenes <- strsplit(as.vector(terms$genes), ', ')
    if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), ',')
    count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
    logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, genes$ID)])
    if(class(logFC) == 'factor'){
      logFC <- gsub(",", ".", gsub("\\.", "", logFC))
      logFC <- as.numeric(logFC)
    }
    s <- 1; zsc <- c()
    for (c in 1:length(count)){
      value <- 0
      e <- s + count[c] - 1
      value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
      zsc <- c(zsc, sum(value) / sqrt(count[c]))
      s <- e + 1
    }
    if (is.null(terms$id)){
      df <- data.frame(category = rep(as.character(terms$category), count), term = rep(as.character(terms$term), count),
                       count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                       zscore = rep(zsc, count), stringsAsFactors = FALSE)
    }else{
      df <- data.frame(category = rep(as.character(terms$category), count), ID = rep(as.character(terms$id), count), term = rep(as.character(terms$term), count),
                       count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                       zscore = rep(zsc, count), stringsAsFactors = FALSE)
    }
    return(df)
  }



  chord_dat <- function(data, genes, process, limit){
    id <- term <- logFC <- BPprocess <- NULL

    if (missing(limit)) limit <- c(0, 0)
    if (missing(genes)){
      if (is.null(data$logFC)){
        genes <- unique(data$genes)
      }else{
        genes <- subset(data, !duplicated(genes), c(genes, logFC))
      }
    }else{
      if(is.vector(genes)){
        genes <- as.character(genes)
      }else{
        if(class(genes[, 2]) != 'numeric') genes[, 2] <- as.numeric(levels(genes[, 2]))[genes[, 2]]
        genes[, 1] <- as.character(genes[, 1])
        colnames(genes) <- c('genes', 'logFC')
      }
    }
    if (missing(process)){
      process <- unique(data$term)
    }else{
      if(class(process) != 'character') process <- as.character(process)
    }
    if (strsplit(process[1],':')[[1]][1] == 'GO'){
      subData <- subset(data, id%in%process)
      colnames(subData)[which(colnames(subData) == 'id')] <- 'BPprocess'
    }else{
      subData <- subset(data, term%in%process)
      colnames(subData)[which(colnames(subData) == 'term')] <- 'BPprocess'
    }

    if(is.vector(genes)){
      M <- genes[genes%in%unique(subData$genes)]
      mat <- matrix(0, ncol = length(process), nrow = length(M))
      rownames(mat) <- M
      colnames(mat) <- process
      for (p in 1:length(process)){
        sub2 <- subset(subData, BPprocess == process[p])
        for (g in 1:length(M)) mat[g, p] <- ifelse(M[g]%in%sub2$genes, 1, 0)
      }
    }else{
      genes <- subset(genes, genes %in% unique(subData$genes))
      N <- length(process) + 1
      M <- genes[,1]
      mat <- matrix(0, ncol = N, nrow = length(M))
      rownames(mat) <- M
      colnames(mat) <- c(process, 'logFC')
      mat[,N] <- genes[,2]
      for (p in 1:(N-1)){
        sub2 <- subset(subData, BPprocess == process[p])
        for (g in 1:length(M)) mat[g, p] <- ifelse(M[g]%in%sub2$genes, 1, 0)
      }
    }
    return(mat)
  }


  GOBubble <- function(data, display, title, color, labels, ID = T, table.legend = T, table.col = T){
    zscore <- adj_pval <- category <- count <- id <- term <- NULL
    if (missing(display)) display <- 'single'
    if (missing(title)) title <- ''
    if (missing(color)) cols <- c("chartreuse4", "brown2", "cornflowerblue") else cols <- color
    if (missing(labels)) labels <- 5

    colnames(data) <- tolower(colnames(data))
    if(!'count'%in%colnames(data)){
      rang <- c(5, 5)
      data$count <- rep(1, dim(data)[1])
    }else {rang <- c(1, 30)}
    data$adj_pval <- -log(data$adj_pval, 10)
    sub <- data[!duplicated(data$term), ]
    g <- ggplot(sub, aes(zscore, adj_pval))+
      labs(title = title, x = 'z-score', y = '-log (adj p-value)')+
      geom_point(aes(col = category, size = count), alpha = 1 / 2)+
      geom_hline(yintercept = 1.3, col = 'orange')+
      scale_size(range = rang, guide = 'none')
    if (!is.character(labels)) sub2 <- subset(sub, subset = sub$adj_pval >= labels) else sub2 <- subset(sub, sub$id%in%labels | sub$term%in%labels)
    if (display == 'single'){
      g <- g + scale_colour_manual('Category', values = cols, labels = c('Biological Process', 'Cellular Component', 'Molecular Function'))+
        theme(legend.position = 'bottom')+
        annotate ("text", x = min(sub$zscore), y = 1.5, label = "threshold", colour = "orange", size = 3)
      if (ID) g <- g+ geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = id), size = 5) else g <- g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = term), size = 4)
      if (table.legend){
        if (table.col) table <- draw_table(sub2, col = cols) else table <- draw_table(sub2)
        g <- g + theme(axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'), panel.background = element_blank(),
                       panel.grid.minor = element_blank(), panel.grid.major = element_line(color = 'grey80'), plot.background = element_blank())
        graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
        grid.arrange(g, table, ncol = 2)
      }else{
        g + theme(axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'), panel.background = element_blank(),
                  panel.grid.minor = element_blank(), panel.grid.major = element_line(color = 'grey80'), plot.background = element_blank())
      }
    }else{
      g <- g + facet_grid(.~category, space = 'free_x', scales = 'free_x') + scale_colour_manual(values = cols, guide ='none')
      if (ID) {
        g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = id), size = 5) +
          theme(axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'), panel.border = element_rect(fill = 'transparent', color = 'grey80'),
                panel.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank())
      }else{
        g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = term), size = 5) +
          theme(axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'), panel.border = element_rect(fill = 'transparent', color = 'grey80'),
                panel.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank())
      }
    }
  }



  GOBar <- function(data, display, order.by.zscore = T, title, zsc.col){
    id <- adj_pval <- zscore <- NULL
    if (missing(display)) display <- 'single'
    if (missing(title)) title <- ''
    if (missing(zsc.col)) zsc.col <- c('red', 'white', 'blue')
    colnames(data) <- tolower(colnames(data))
    data$adj_pval <- -log(data$adj_pval, 10)
    sub <- data[!duplicated(data$term), ]

    if (order.by.zscore == T) {
      sub <- sub[order(sub$zscore, decreasing = T), ]
      leg <- theme(legend.position = 'bottom')
      g <-  ggplot(sub, aes(x = factor(id, levels = stats::reorder(id, adj_pval)), y = adj_pval, fill = zscore)) +
        geom_bar(stat = 'identity', color = 'black') +
        scale_fill_gradient2('z-score', low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
                             breaks = c(min(sub$zscore), max(sub$zscore)), labels = c('decreasing', 'increasing')) +
        labs(title = title, x = '', y = '-log (adj p-value)') +
        leg
    }else{
      sub <- sub[order(sub$adj_pval, decreasing = T), ]
      leg <- theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.995), legend.background = element_rect(fill = 'transparent'),
                   legend.box = 'vertical', legend.direction = 'horizontal')
      g <-  ggplot(sub, aes( x = factor(id, levels = reorder(id, adj_pval)), y = zscore, fill = adj_pval)) +
        geom_bar(stat = 'identity', color = 'black') +
        scale_fill_gradient2('Significance', guide = guide_colorbar(title.position = "top", title.hjust = 0.5), breaks = c(min(sub$adj_pval), max(sub$adj_pval)), labels = c('low', 'high')) +
        labs(title = title, x = '', y = 'z-score') +
        leg
    }
    if (display == 'single'){
      g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'),
                panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), plot.background = element_blank())
    }else{
      g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.line = element_line(color = 'grey80'), axis.ticks = element_line(color = 'grey80'),
                panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), plot.background = element_blank())+
        facet_grid(.~category, space = 'free_x', scales = 'free_x')
    }
  }



  GOCircle <- function(data, title, nsub, rad1, rad2, table.legend = T, zsc.col, lfc.col, label.size, label.fontface){
    xmax <- y1<- zscore <- y2 <- ID <- logx <- logy2 <- logy <- logFC <- NULL
    if (missing(title)) title <- ''
    if (missing(nsub)) if (dim(data)[1] > 10) nsub <- 10 else nsub <- dim(data)[1]
    if (missing(rad1)) rad1 <- 2
    if (missing(rad2)) rad2 <- 3
    if (missing(zsc.col)) zsc.col <- c('red', 'white', 'blue')
    if (missing(lfc.col)) lfc.col <- c('cornflowerblue', 'firebrick1') else lfc.col <- rev(lfc.col)
    if (missing(label.size)) label.size = 5
    if (missing(label.fontface)) label.fontface = 'bold'

    data$adj_pval <- -log(data$adj_pval, 10)
    suby <- data[!duplicated(data$term), ]
    if (is.numeric(nsub) == T){
      suby <- suby[1:nsub, ]
    }else{
      if (strsplit(nsub[1], ':')[[1]][1] == 'GO'){
        suby <- suby[suby$ID%in%nsub, ]
      }else{
        suby <- suby[suby$term%in%nsub, ]
      }
      nsub <- length(nsub)}
    N <- dim(suby)[1]
    r_pval <- round(range(suby$adj_pval), 0) + c(-2, 2)
    ymax <- c()
    for (i in 1:length(suby$adj_pval)){
      val <- (suby$adj_pval[i] - r_pval[1]) / (r_pval[2] - r_pval[1])
      ymax <- c(ymax, val)}
    df <- data.frame(x = seq(0, 10 - (10 / N), length = N), xmax = rep(10 / N - 0.2, N), y1 = rep(rad1, N), y2 = rep(rad2, N), ymax = ymax, zscore = suby$zscore, ID = suby$ID)
    scount <- data[!duplicated(data$term), which(colnames(data) == 'count')][1:nsub]
    idx_term <- which(!duplicated(data$term) == T)
    xm <- c(); logs <- c()
    for (sc in 1:length(scount)){
      idx <- c(idx_term[sc], idx_term[sc + 1] - 1)
      val <- stats::runif(scount[sc], df$x[sc] + 0.06, (df$x[sc] + df$xmax[sc] - 0.06))
      xm <- c(xm, val)
      r_logFC <- round(range(data$logFC[idx[1]:idx[2]]), 0) + c(-1, 1)
      for (lfc in idx[1]:idx[2]){
        val <- (data$logFC[lfc] - r_logFC[1]) / (r_logFC[2] - r_logFC[1])
        logs <- c(logs, val)}
    }
    cols <- c()
    for (ys in 1:length(logs)) cols <- c(cols, ifelse(data$logFC[ys] > 0, 'upregulated', 'downregulated'))
    dfp <- data.frame(logx = xm, logy = logs, logFC = factor(cols), logy2 = rep(rad2, length(logs)))
    c <-  ggplot()+
      geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y1, ymax = y1 + ymax, fill = zscore), colour = 'black') +
      geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y2, ymax = y2 + 1), fill = 'gray70') +
      geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y2 + 0.5, ymax = y2 + 0.5), colour = 'white') +
      geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y2 + 0.25, ymax = y2 + 0.25), colour = 'white') +
      geom_rect(data = df, aes(xmin = x, xmax = x + xmax, ymin = y2 + 0.75, ymax = y2 + 0.75), colour = 'white') +
      geom_text(data = df, aes(x = x + (xmax / 2), y = y2 + 1.3, label = ID, angle = 360 - (x = x + (xmax / 2)) / (10 / 360)), size = label.size, fontface = label.fontface) +
      coord_polar() +
      labs(title = title) +
      ylim(1, rad2 + 1.6) +
      xlim(0, 10) +
      theme_blank +
      scale_fill_gradient2('z-score', low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], guide = guide_colorbar(title.position = "top", title.hjust = 0.5), breaks = c(min(df$zscore), max(df$zscore)),labels = c('decreasing', 'increasing')) +
      theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal') +
      geom_point(data = dfp, aes(x = logx, y = logy2 + logy), pch = 21, fill = 'transparent', color = 'black', size = 3)+
      geom_point(data = dfp, aes(x = logx, y = logy2 + logy, color = logFC), size = 2.5)+
      scale_colour_manual(values = lfc.col, guide = guide_legend(title.position = "top", title.hjust = 0.5))

    if (table.legend){
      table <- draw_table(suby, col = 'black')
      graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
      grid.arrange(c, table, ncol = 2)
    }else{
      c + theme(plot.background = element_rect(fill = 'aliceblue'), panel.background = element_rect(fill = 'white'))
    }
  }



  #source("C:/Users/Yang/Documents/mmease001/R/GOVenn.R")


  circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }


  get_overlap<-function(A,B,C){
    colnames(A)<-c('ID','logFC')
    colnames(B)<-c('ID','logFC')
    colnames(C)<-c('ID','logFC')
    UP<-NULL;DOWN<-NULL;Change<-NULL
    if (class(A$logFC)!='numeric'){
      A$logFC<-gsub(",", ".", gsub("\\.", "", A$logFC))
      A$Trend<-sapply(as.numeric(A$logFC), function(x) ifelse(x > 0,'UP','DOWN'))
    }else{ A$Trend<-sapply(A$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
    if (class(B$logFC)!='numeric'){
      B$logFC<-gsub(",", ".", gsub("\\.", "", B$logFC))
      B$Trend<-sapply(as.numeric(B$logFC), function(x) ifelse(x > 0,'UP','DOWN'))
    }else{ B$Trend<-sapply(B$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
    if (class(C$logFC)!='numeric'){
      C$logFC<-gsub(",", ".", gsub("\\.", "", C$logFC))
      C$Trend<-sapply(as.numeric(C$logFC), function(x) ifelse(x > 0,'UP','DOWN'))
    }else{ C$Trend<-sapply(C$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
    if (sum(((A$ID%in%B$ID)==T)==T)==0){
      AB<-data.frame()
    }else{
      AB<-A[(A$ID%in%B$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
      BA<-B[(B$ID%in%A$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
      AB<-merge(AB,BA,by="ID")
      rownames(AB)<-AB$ID
      AB<-AB[,-1]
    }
    if (sum(((A$ID%in%C$ID)==T)==T)==0){
      AC<-data.frame()
    }else{
      AC<-A[(A$ID%in%C$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
      CA<-C[(C$ID%in%A$ID)==T,which(colnames(C)%in%c('ID','logFC','Trend'))]
      AC<-merge(AC,CA,by="ID")
      rownames(AC)<-AC$ID
      AC<-AC[,-1]
    }
    if (sum(((B$ID%in%C$ID)==T)==T)==0){
      BC<-data.frame()
    }else{
      BC<-B[(B$ID%in%C$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
      CB<-C[(C$ID%in%B$ID)==T,which(colnames(C)%in%c('ID','logFC','Trend'))]
      BC<-merge(BC,CB,by="ID")
      rownames(BC)<-BC$ID
      BC<-BC[,-1]
    }
    if (sum(((A$ID%in%B$ID)==T & (A$ID%in%C$ID)==T))==0){
      ABC<-data.frame()
    }else{
      ABC<-A[((A$ID%in%B$ID)==T & (A$ID%in%C$ID)==T),which(colnames(A)%in%c('ID','logFC','Trend'))]
      BAC<-B[((B$ID%in%A$ID)==T & (B$ID%in%C$ID)==T),which(colnames(B)%in%c('ID','logFC','Trend'))]
      CAB<-C[((C$ID%in%A$ID)==T & (C$ID%in%B$ID)==T),which(colnames(C)%in%c('ID','logFC','Trend'))]
      ABC<-merge(ABC,BAC,by='ID')
      ABC<-merge(ABC,CAB,by='ID')
      rownames(ABC)<-ABC$ID
      ABC<-ABC[,-1]
    }
    A_only<-A[((A$ID%in%B$ID)==F & (A$ID%in%C$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
    rownames(A_only)<-A_only$ID
    A_only<-A_only[,-1]
    B_only<-B[((B$ID%in%A$ID)==F & (B$ID%in%C$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
    rownames(B_only)<-B_only$ID
    B_only<-B_only[,-1]
    C_only<-C[((C$ID%in%A$ID)==F & (C$ID%in%B$ID)==F),which(colnames(A)%in%c('ID','logFC','Trend'))]
    rownames(C_only)<-C_only$ID
    C_only<-C_only[,-1]
    UP<-c(UP,sum(A_only$Trend=='UP'));DOWN<-c(DOWN,sum(A_only$Trend=='DOWN'));Change<-c(Change,sum(A_only$Trend=='Change'))
    UP<-c(UP,sum(B_only$Trend=='UP'));DOWN<-c(DOWN,sum(B_only$Trend=='DOWN'));Change<-c(Change,sum(B_only$Trend=='Change'))
    UP<-c(UP,sum(C_only$Trend=='UP'));DOWN<-c(DOWN,sum(C_only$Trend=='DOWN'));Change<-c(Change,sum(C_only$Trend=='Change'))
    if (dim(AB)[1]==0){
      OvAB<-data.frame()
      UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
    }else{
      tmp<-NULL
      for (t in 1:dim(AB)[1]) tmp<-c(tmp,ifelse(AB$Trend.x[t]==AB$Trend.y[t],AB$Trend.x[t],'Change'))
      OvAB<-data.frame(logFC_A=AB$logFC.x,logFC_B=AB$logFC.y,Trend=tmp)
      rownames(OvAB)<-rownames(AB)
      AB<-OvAB[order(OvAB$Trend),]
      UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
    }
    if (dim(AC)[1]==0){
      OvAc<-data.frame()
      UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
    }else{
      tmp<-NULL
      for (t in 1:dim(AC)[1]) tmp<-c(tmp,ifelse(AC$Trend.x[t]==AC$Trend.y[t],AC$Trend.x[t],'Change'))
      OvAC<-data.frame(logFC_A=AC$logFC.x,logFC_C=AC$logFC.y,Trend=tmp)
      rownames(OvAC)<-rownames(AC)
      AC<-OvAC[order(OvAC$Trend),]
      UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
    }
    if (dim(BC)[1]==0){
      OvBC<-data.frame()
      UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
    }else{
      tmp<-NULL
      for (t in 1:dim(BC)[1]) tmp<-c(tmp,ifelse(BC$Trend.x[t]==BC$Trend.y[t],BC$Trend.x[t],'Change'))
      OvBC<-data.frame(logFC_B=BC$logFC.x,logFC_C=BC$logFC.y,Trend=tmp)
      rownames(OvBC)<-rownames(BC)
      BC<-OvBC[order(OvBC$Trend),]
      UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
    }
    if (dim(ABC)[1]==0){
      OvABC<-data.frame()
      UP<-c(UP,0);DOWN<-c(DOWN,0);Change<-c(Change,0)
    }else{
      tmp<-NULL
      for (t in 1:dim(ABC)[1]) tmp<-c(tmp,ifelse(((ABC$Trend.x[t]==ABC$Trend.y[t]) & (ABC$Trend.x[t]==ABC$Trend[t])),ABC$Trend.x[t],'Change'))
      OvABC<-data.frame(logFC_A=ABC$logFC.x,logFC_B=ABC$logFC.y,logFC_C=ABC$logFC,Trend=tmp)
      rownames(OvABC)<-rownames(ABC)
      ABC<-OvABC[order(OvABC$Trend),]
      UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
    }
    counts<-data.frame(Contrast=c('A_only','B_only','C_only','AB','AC','BC','ABC'),Count=c(dim(A_only)[1],dim(B_only)[1],dim(C_only)[1],dim(AB)[1],dim(AC)[1],dim(BC)[1],dim(ABC)[1]),UP=UP,DOWN=DOWN,Change=Change)
    venn<-list(A_only=A_only,B_only=B_only,C_only=C_only,AB=AB,BC=BC,AC=AC,ABC=ABC)
    return(list(venn_df=counts,table=venn))
  }


  get_overlap2<-function(A,B){
    colnames(A)<-c('ID','logFC')
    colnames(B)<-c('ID','logFC')
    UP<-NULL;DOWN<-NULL;Change<-NULL
    if (class(A$logFC)!='numeric'){
      A$logFC<-gsub(",", ".", gsub("\\.", "", A$logFC))
      A$Trend<-sapply(as.numeric(A$logFC), function(x) ifelse(x > 0,'UP','DOWN'))
    }else{ A$Trend<-sapply(A$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
    if (class(B$logFC)!='numeric'){
      B$logFC<-gsub(",", ".", gsub("\\.", "", B$logFC))
      B$Trend<-sapply(as.numeric(B$logFC), function(x) ifelse(x > 0,'UP','DOWN'))
    }else{ B$Trend<-sapply(B$logFC, function(x) ifelse(x > 0,'UP','DOWN'))}
    AB<-A[(A$ID%in%B$ID)==T,which(colnames(A)%in%c('ID','logFC','Trend'))]
    BA<-B[(B$ID%in%A$ID)==T,which(colnames(B)%in%c('ID','logFC','Trend'))]
    A_only<-A[(A$ID%in%B$ID)==F,which(colnames(A)%in%c('ID','logFC','Trend'))]
    B_only<-B[(B$ID%in%A$ID)==F,which(colnames(B)%in%c('ID','logFC','Trend'))]
    AB<-merge(AB,BA,by='ID')
    UP<-c(UP,sum(A_only$Trend=='UP'));DOWN<-c(DOWN,sum(A_only$Trend=='DOWN'));Change<-c(Change,sum(A_only$Trend=='Change'))
    UP<-c(UP,sum(B_only$Trend=='UP'));DOWN<-c(DOWN,sum(B_only$Trend=='DOWN'));Change<-c(Change,sum(B_only$Trend=='Change'))
    rownames(A_only)<-A_only$ID
    A_only<-A_only[,-1]
    A_only<-A_only[order(A_only$Trend),]
    rownames(B_only)<-B_only$ID
    B_only<-B_only[,-1]
    B_only<-B_only[order(B_only$Trend),]
    tmp<-NULL
    for (t in 1:dim(AB)[1]) tmp<-c(tmp,ifelse(AB$Trend.x[t]==AB$Trend.y[t],AB$Trend.x[t],'Change'))
    OvAB<-data.frame(logFC_A=AB$logFC.x,logFC_B=AB$logFC.y,Trend=tmp)
    rownames(OvAB)<-AB$ID
    AB<-OvAB[order(OvAB$Trend),]
    UP<-c(UP,sum(tmp=='UP'));DOWN<-c(DOWN,sum(tmp=='DOWN'));Change<-c(Change,sum(tmp=='Change'))
    counts<-data.frame(Contrast=c('A_only','B_only','AB'),Count=c(dim(A_only)[1],dim(B_only)[1],dim(AB)[1]),UP=UP,DOWN=DOWN,Change=Change)
    venn<-list(A_only=A_only,B_only=B_only,AB=AB)
    return(list(venn_df=counts,table=venn,dim=c(dim(A)[1],dim(B)[1])))
  }



  GOVenn<-function(data1, data2, data3, title, label, lfc.col, circle.col, plot=T){
    id <- NULL
    if (missing(label)) label<-c('List1','List2','List3')
    if (missing(lfc.col)) lfc.col<-c('firebrick1','gold','cornflowerblue')
    if (missing(circle.col)) circle.col<-c('brown1','chartreuse3','cornflowerblue')
    if (missing(title)) title<-''
    if (missing(data3)==F) {
      three<-T
      overlap<-get_overlap(data1,data2,data3)
      venn_df<-overlap$venn_df
      table<-overlap$table
    }else{
      three<-F
      overlap<-get_overlap2(data1,data2)
      venn_df<-overlap$venn_df
      table<-overlap$table
    }

    ### calc Venn ###
    if (three){
      center<-data.frame(x=c(0.4311,0.4308,0.6380),y=c(0.6197,0.3801,0.5001),diameter=c(0.4483,0.4483,0.4483))
      outerCircle<-data.frame(x=numeric(),y=numeric(),id=numeric())
      for (var in 1:3){
        dat <- circleFun(c(center$x[var],center$y[var]),center$diameter[var],npoints = 100)
        outerCircle<-rbind(outerCircle,dat)
      }
      outerCircle$id<-rep(c(label[1],label[2],label[3]),each=100)
      outerCircle$id<-factor(outerCircle$id, levels=c(label[1],label[2],label[3]))
    }else{
      center<-data.frame(x=c(0.33,0.6699),y=c(0.5,0.5),diameter=c(0.6180,0.6180))
      outerCircle<-data.frame(x=numeric(),y=numeric(),id=numeric())
      for (var in 1:2){
        dat <- circleFun(c(center$x[var],center$y[var]),center$diameter[var],npoints = 100)
        outerCircle<-rbind(outerCircle,dat)
      }
      outerCircle$id<-rep(c(label[1],label[2]),each=100)
      outerCircle$id<-factor(outerCircle$id, levels=c(label[1],label[2]))
    }

    ### calc single pies ###
    if (three){
      Pie<-data.frame(x=numeric(),y=numeric(),id=numeric())
      dat <- circleFun(c(center$x[1],max(subset(outerCircle,id==label[1])$y)-0.05),0.1,npoints = 100)
      Pie<-rbind(Pie,dat)
      dat <- circleFun(c(center$x[2],min(subset(outerCircle,id==label[2])$y)+0.05),0.1,npoints = 100)
      Pie<-rbind(Pie,dat)
      dat <- circleFun(c(max(subset(outerCircle,id==label[3])$x)-0.05,center$y[3]),0.1,npoints = 100)
      Pie<-rbind(Pie,dat)
      Pie$id<-rep(1:3,each=100)
      UP<-Pie[c(1:50,100:150,200:250),]
      Down<-Pie[c(50:100,150:200,250:300),]
    }else{
      Pie<-data.frame(x=numeric(),y=numeric(),id=numeric())
      dat <- circleFun(c(min(subset(outerCircle,id==label[1])$x)+0.05,center$y[1]),0.1,npoints = 100)
      Pie<-rbind(Pie,dat)
      dat <- circleFun(c(max(subset(outerCircle,id==label[2])$x)-0.05,center$y[2]),0.1,npoints = 100)
      Pie<-rbind(Pie,dat)
      Pie$id<-rep(1:2,each=100)
      UP<-Pie[c(1:50,100:150),]
      Down<-Pie[c(50:100,150:200),]
    }

    ### calc single pie text ###
    if (three){
      x<-c();y<-c()
      for (i in unique(Pie$id)){
        x<-c(x,rep((min(subset(Pie,id==i)$x)+max(subset(Pie,id==i)$x))/2,2))
        y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2+0.02)
        y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2-0.02)
      }
      pieText<-data.frame(x=x,y=y,label=c(venn_df$UP[1],venn_df$DOWN[1],venn_df$UP[2],venn_df$DOWN[2],venn_df$UP[3],venn_df$DOWN[3]))
    }else{
      x<-c();y<-c()
      for (i in unique(Pie$id)){
        x<-c(x,rep((min(subset(Pie,id==i)$x)+max(subset(Pie,id==i)$x))/2,2))
        y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2+0.02)
        y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2-0.02)
      }
      pieText<-data.frame(x=x,y=y,label=c(venn_df$UP[1],venn_df$DOWN[1],venn_df$UP[2],venn_df$DOWN[2]))
    }

    ### calc overlap pies ###
    if (three){
      smc<-data.frame(x=c(0.6,0.59,0.31,0.5),y=c(0.66,0.34,0.5,0.5))
      PieOv<-data.frame(x=numeric(),y=numeric())
      PieOv<-rbind(PieOv,circleFun(c(smc$x[1],smc$y[1]),0.06,npoints = 100))
      PieOv<-rbind(PieOv,circleFun(c(smc$x[2],smc$y[2]),0.06,npoints = 100))
      PieOv<-rbind(PieOv,circleFun(c(smc$x[3],smc$y[3]),0.06,npoints = 100))
      PieOv<-rbind(PieOv,circleFun(c(smc$x[4],smc$y[4]),0.06,npoints = 100))
      PieOv$id<-rep(1:4,each=100)
      smc$id<-1:4
      UPOv<-rbind(smc[1,],PieOv[1:33,],smc[1,],smc[2,],PieOv[100:133,],smc[2,],smc[3,],PieOv[200:233,],smc[3,],smc[4,],PieOv[300:333,],smc[4,])
      Change<-rbind(smc[1,],PieOv[33:66,],smc[1,],smc[2,],PieOv[133:166,],smc[2,],smc[3,],PieOv[233:266,],smc[3,],smc[4,],PieOv[333:366,],smc[4,])
      DownOv<-rbind(smc[1,],PieOv[66:100,],smc[1,],smc[2,],PieOv[166:200,],smc[2,],smc[3,],PieOv[266:300,],smc[3,],smc[4,],PieOv[366:400,],smc[4,])
    }else{
      PieOv<-data.frame(x=numeric(),y=numeric(),id=numeric())
      PieOv<-rbind(PieOv,circleFun(c(0.5,0.5),0.08,npoints = 100))
      PieOv$id<-rep(1,100)
      center<-data.frame(x=0.5, y=0.5, id=1)
      UPOv<-rbind(center[1,],PieOv[1:33,])
      Change<-rbind(center[1,],PieOv[33:66,])
      DownOv<-rbind(center[1,],PieOv[66:100,])
    }

    ### calc overlap pie text ###
    if (three){
      x<-c();y<-c()
      for (i in unique(PieOv$id)){
        x<-c(x,subset(UPOv,id==i)$x[1]+0.0115,subset(DownOv,id==i)$x[1]-0.018,subset(Change,id==i)$x[1]+0.01)
        y<-c(y,subset(UPOv,id==i)$y[1]+0.01,subset(DownOv,id==i)$y[1],subset(Change,id==i)$y[1]-0.013)
      }
      small.pieT<-data.frame(x=x,y=y,label=c(venn_df$UP[5],venn_df$Change[5],venn_df$DOWN[5],venn_df$UP[6],venn_df$Change[6],venn_df$DOWN[6],venn_df$UP[4],venn_df$Change[4],venn_df$DOWN[4],venn_df$UP[7],venn_df$Change[7],venn_df$DOWN[7]))
    }else{
      x<-c(subset(UPOv,id==1)$x[1]+0.015,subset(DownOv,id==1)$x[1]-0.018,subset(Change,id==1)$x[1]+0.01)
      y<-c(subset(UPOv,id==1)$y[1]+0.015,subset(DownOv,id==1)$y[1],subset(Change,id==1)$y[1]-0.013)
      small.pieT<-data.frame(x=x,y=y,label=c(venn_df$UP[3],venn_df$Change[3],venn_df$DOWN[3]))
    }

    g<- ggplot()+
      geom_polygon(data=outerCircle, aes(x,y, group=id, fill=id) ,alpha=0.5,color='black')+
      scale_fill_manual(values=circle.col)+
      guides(fill=guide_legend(title=''))+
      geom_polygon(data=UP, aes(x,y,group=id),fill=lfc.col[1],color='white')+
      geom_polygon(data=Down, aes(x,y,group=id),fill=lfc.col[3],color='white')+
      geom_text(data=pieText, aes(x=x,y=y,label=label),size=5)+
      geom_polygon(data=UPOv, aes(x,y,group=id),fill=lfc.col[1],color='white')+
      geom_polygon(data=DownOv, aes(x,y,group=id),fill=lfc.col[3],color='white')+
      geom_polygon(data=Change, aes(x,y,group=id),fill=lfc.col[2],color='white')+
      geom_text(data=small.pieT,aes(x=x,y=y,label=label),size=4)+
      theme_blank+
      labs(title=title)

    if (plot) return(g) else return(list(plot=g,table=table))
  }



  #source("C:/Users/Yang/Documents/mmease001/R/GOCluster.R")


  GOCluster<-function(data, process, metric, clust, clust.by, nlfc, lfc.col, lfc.min, lfc.max, lfc.space, lfc.width, term.col, term.space, term.width){
    x <- y <- xend <- yend <- width <- space <- logFC <- NULL
    if (missing(metric)) metric<-'euclidean'
    if (missing(clust)) clust<-'average'
    if (missing(clust.by)) clust.by<-'term'
    if (missing(nlfc)) nlfc<-F
    if (missing(lfc.col)) lfc.col<-c('firebrick1','white','dodgerblue')
    if (missing(lfc.min)) lfc.min <- -3
    if (missing(lfc.max)) lfc.max <- 3
    if (missing(lfc.space)) lfc.space<- (-0.5) else lfc.space<-lfc.space*(-1)
    if (missing(lfc.width)) lfc.width<- (-1.6) else lfc.width<-lfc.space-lfc.width-0.1
    if (missing(term.col)) term.col<-brewer.pal(length(process), 'Set3')
    if (missing(term.space)) term.space<- lfc.space+lfc.width else term.space<-term.space*(-1)+lfc.width
    if (missing(term.width)) term.width<- 2*lfc.width+term.space else term.width<-term.width*(-1)+term.space

    if (nlfc){
      colnames(data)[1:3] <- c('genes','term','logFC')
      chord <- chord_dat(data[,1:3])
    }else{
      chord <- chord_dat(data = data, process = process)
    }
    if (clust.by=='logFC') distance <- stats::dist(chord[,dim(chord)[2]], method=metric)
    if (clust.by=='term') distance <- stats::dist(chord, method=metric)
    cluster <- stats::hclust(distance, method=clust)
    dendr <- dendro_data(cluster)
    y_range <- range(dendr$segments$y)
    x_pos <- data.frame(x=dendr$label$x, label=as.character(dendr$label$label))
    chord <- as.data.frame(chord)
    chord$label <- as.character(rownames(chord))
    all <- merge(x_pos, chord, by='label')
    all$label <- as.character(all$label)
    if (nlfc){
      lfc_rect <- all[,c(2, dim(all)[2])]
      for (l in 4:dim(data)[2]) lfc_rect <- cbind(lfc_rect, sapply(all$label, function(x) data[match(x, data$genes), l]))
      num <- dim(data)[2]-1
      tmp <- seq(lfc.space, lfc.width, length = num)
      lfc<-data.frame(x=numeric(),width=numeric(),space=numeric(),logFC=numeric())
      for (l in 1:(length(tmp)-1)){
        tmp_df<-data.frame(x=lfc_rect[,1],width=tmp[l+1],space=tmp[l],logFC=lfc_rect[,l+1])
        lfc<-rbind(lfc,tmp_df)
      }
    }else{
      lfc <- all[,c(2, dim(all)[2])]
      lfc$space <- lfc.space
      lfc$width <- lfc.width
    }
    term <- all[,c(2:(length(process)+2))]
    color<-NULL;termx<-NULL;tspace<-NULL;twidth<-NULL
    for (row in 1:dim(term)[1]){
      idx <- which(term[row,-1] != 0)
      if(length(idx) != 0){
        termx<-c(termx,rep(term[row,1],length(idx)))
        color<-c(color,term.col[idx])
        tmp<-seq(term.space,term.width,length=length(idx)+1)
        tspace<-c(tspace,tmp[1:(length(tmp)-1)])
        twidth<-c(twidth,tmp[2:length(tmp)])
      }
    }
    tmp <- sapply(lfc$logFC, function(x) ifelse(x > lfc.max, lfc.max, x))
    logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, lfc.min, x))
    lfc$logFC <- logFC
    term_rect <- data.frame(x = termx, width = twidth, space = tspace, col = color)
    legend <- data.frame(x = 1:length(process),label = process)

    ggplot()+
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend))+
      geom_rect(data=lfc,aes(xmin=x-0.5,xmax=x+0.5,ymin=width,ymax=space,fill=logFC))+
      scale_fill_gradient2('logFC',low=lfc.col[3],mid=lfc.col[2],high=lfc.col[1],guide=guide_colorbar(title.position='top',title.hjust=0.5),breaks=c(min(lfc$logFC),max(lfc$logFC)),labels=c(round(min(lfc$logFC)),round(max(lfc$logFC))))+
      geom_rect(data=term_rect,aes(xmin=x-0.5,xmax=x+0.5,ymin=width,ymax=space),fill=term_rect$col)+
      geom_point(data=legend,aes(x=x,y=0.1,size=factor(label,levels=label),shape=NA))+
      guides(size=guide_legend("Terms",ncol=4,byrow=T,override.aes=list(shape=22,fill=term.col,size = 8)))+
      coord_polar()+
      scale_y_reverse()+
      theme(legend.position='bottom',legend.background = element_rect(fill='transparent'),legend.box='horizontal',legend.direction='horizontal')+
      theme_blank

  }



  GOChord <- function(data, title, space, gene.order, gene.size, gene.space, nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size, process.label){
    y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
    Ncol <- dim(data)[2]

    if (missing(title)) title <- ''
    if (missing(space)) space = 0
    if (missing(gene.order)) gene.order <- 'none'
    if (missing(gene.size)) gene.size <- 3
    if (missing(gene.space)) gene.space <- 0.2
    if (missing(lfc.col)) lfc.col <- c('brown1', 'azure', 'cornflowerblue')
    if (missing(lfc.min)) lfc.min <- -3
    if (missing(lfc.max)) lfc.max <- 3
    if (missing(ribbon.col)) colRib <- grDevices::rainbow(Ncol - nlfc) else colRib <- ribbon.col
    if (missing(border.size)) border.size <- 0.5
    if (missing (process.label)) process.label <- 11

    if (gene.order == 'logFC') data <- data[order(data[, Ncol], decreasing = T), ]
    if (gene.order == 'alphabetical') data <- data[order(rownames(data)), ]
    if (sum(!is.na(match(colnames(data), 'logFC'))) > 0){
      if (nlfc == 1){
        cdata <- data[, 1:(Ncol - 1)]
        lfc <-data[, Ncol]
      }else{
        cdata <- data[, 1:(Ncol - nlfc)]
        lfc <- data[, (Ncol - nlfc + 1)]
      }
    }else{
      cdata <- data
      lfc <- 0
    }
    nrib <- colSums(cdata)
    ngen <- rowSums(cdata)
    Ncol <- dim(cdata)[2]
    Nrow <- dim(cdata)[1]
    colRibb <- c()
    for (b in 1:length(nrib)) colRibb <- c(colRibb, rep(colRib[b], 202 * nrib[b]))
    r1 <- 1; r2 <- r1 + 0.1
    xmax <- c(); x <- 0
    for (r in 1:length(nrib)){
      perc <- nrib[r] / sum(nrib)
      xmax <- c(xmax, (pi * perc) - space)
      if (length(x) <= Ncol - 1) x <- c(x, x[r] + pi * perc)
    }
    xp <- c(); yp <- c()
    l <- 50
    for (s in 1:Ncol){
      xh <- seq(x[s], x[s] + xmax[s], length = l)
      xp <- c(xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] + xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)), r2 * sin(x[s]))
      yp <- c(yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] + xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)), r2 * cos(x[s]))
    }
    df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol), each = 4 + 2 * l))
    xp <- c(); yp <- c(); logs <- NULL
    x2 <- seq(0 - space, -pi - (-pi / Nrow) - space, length = Nrow)
    xmax2 <- rep(-pi / Nrow + space, length = Nrow)
    for (s in 1:Nrow){
      xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
      if (nlfc <= 1){
        xp <- c(xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) * sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]), r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)), r2 * sin(x2[s]))
        yp <- c(yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) * cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]), r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)), r2 * cos(x2[s]))
      }else{
        tmp <- seq(r1, r2, length = nlfc + 1)
        for (t in 1:nlfc){
          logs <- c(logs, data[s, (dim(data)[2] + 1 - t)])
          xp <- c(xp, (tmp[t]) * sin(x2[s]), (tmp[t]) * sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]), tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t + 1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s]))
          yp <- c(yp, (tmp[t]) * cos(x2[s]), (tmp[t]) * cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]), tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t + 1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s]))
        }}}
    if(lfc[1] != 0){
      if (nlfc == 1){
        df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), each = 4 + 2 * l), logFC = rep(lfc, each = 4 + 2 * l))
      }else{
        df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc*Nrow)), each = 4 + 2 * l), logFC = rep(logs, each = 4 + 2 * l))
      }
    }else{
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), each = 4 + 2 * l))
    }
    aseq <- seq(0, 180, length = length(x2)); angle <- c()
    for (o in aseq) if((o + 270) <= 360) angle <- c(angle, o + 270) else angle <- c(angle, o - 90)
    df_texg <- data.frame(xgen = (r1 + gene.space) * sin(x2 + xmax2/2),ygen = (r1 + gene.space) * cos(x2 + xmax2 / 2),labels = rownames(cdata), angle = angle)
    df_texp <- data.frame(xpro = (r1 + 0.15) * sin(x + xmax / 2),ypro = (r1 + 0.15) * cos(x + xmax / 2), labels = colnames(cdata), stringsAsFactors = FALSE)
    cols <- rep(colRib, each = 4 + 2 * l)
    x.end <- c(); y.end <- c(); processID <- c()
    for (gs in 1:length(x2)){
      val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] + 1)
      pros <- which((cdata[gs, ] != 0) == T)
      for (v in 1:(length(val) - 1)){
        x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
        y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
        processID <- c(processID, rep(pros[v], 2))
      }
    }
    df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
    df_bezier <- df_bezier[order(df_bezier$processID,-df_bezier$y.end),]
    x.start <- c(); y.start <- c()
    for (rs in 1:length(x)){
      val<-seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] + 1)
      for (v in 1:(length(val) - 1)){
        x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
        y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
      }
    }
    df_bezier$x.start <- x.start
    df_bezier$y.start <- y.start
    df_path <- bezier(df_bezier, colRib)
    if(length(df_genes$logFC) != 0){
      tmp <- sapply(df_genes$logFC, function(x) ifelse(x > lfc.max, lfc.max, x))
      logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, lfc.min, x))
      df_genes$logFC <- logFC
    }

    g<- ggplot() +
      geom_polygon(data = df_process, aes(x, y, group=id), fill='gray70', inherit.aes = F,color='black') +
      geom_polygon(data = df_process, aes(x, y, group=id), fill=cols, inherit.aes = F,alpha=0.6,color='black') +
      geom_point(aes(x = xpro, y = ypro, size = factor(labels, levels = labels), shape = NA), data = df_texp) +
      #guides(size = guide_legend("Enriched Terms", ncol = 4, byrow = T, override.aes = list(shape = 22, fill = unique(cols), size = 8))) +
      guides(size = guide_legend("", ncol = 3, byrow = T, override.aes = list(shape = 22, fill = unique(cols), size = 8))) +
      theme(legend.text = element_text(size = process.label)) +
      geom_text(aes(xgen, ygen, label = labels, angle = angle), data = df_texg, size = gene.size) +
      geom_polygon(aes(x = lx, y = ly, group = ID), data = df_path, fill = colRibb, color = 'black', size = border.size, inherit.aes = F) +
      labs(title = title) +
      theme_blank

    if (nlfc >= 1){
      g + geom_polygon(data = df_genes, aes(x, y, group = id, fill = logFC), inherit.aes = F, color = 'black') +
        scale_fill_gradient2('logFC', low = lfc.col[3], mid = lfc.col[2], high = lfc.col[1], guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
                             breaks = c(min(df_genes$logFC), max(df_genes$logFC)), labels = c(round(min(df_genes$logFC)), round(max(df_genes$logFC)))) +
        theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal')
    }else{
      g + geom_polygon(data = df_genes, aes(x, y, group = id), fill = 'gray50', inherit.aes = F, color = 'black')+
        theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal')
    }
  }




  GOChord(dat.chord, lfc.col=c('gray0','gray25','gray50'), space=0.02, process.label=11)
  #GOChord(dat.chord, lfc.col=c('gray0','gray25','gray50'), lfc.min=input$KEGGPTrange[1], lfc.max=input$KEGGPTrange[2])

}



plot_bar <- function(matrix){


  pv_col <- function(x){
    1-(x-min(x))/(max(x)-min(x))
  }

  pv <- pv_col(-log10(matrix["P.Val"]))

  cl <- c()
  for (i in 1:length(matrix$P.Val)){
    cl[i]=rgb(1,pv[i,1],0)
  }

  p1 <- ggplot(data=matrix, aes(x=reorder(Category.name, -P.Val), y=Count)) + theme_bw() + coord_flip()

  p2 <- p1 + geom_bar(stat="identity", width = 0.5,aes(fill=factor(P.Val))) +
    scale_fill_manual(values=cl, name="P value")

  p3 <- p2 + theme(#axis.line.x.bottom = element_line(),
    axis.text.x = element_text(color = "black",size=16),
    axis.text.y = element_text(color = "black",size=16),
    axis.title.x = element_text(color = "black",size=16),
    axis.title.y=element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank())

  p4 <- p3 + labs(title = "Enrichment Overview") +
    theme(plot.title = element_text(size = 20, color = "black", hjust = 0, face = "bold"),
          axis.line.x.bottom = element_line())

  p4

}




#' Title
#'
#' @param db The database of input for metabolite enrichment.
#' @param EnrichResultList The table of results for metabolite enrichment.
#'
#' @return The plot of results for metabolite enrichment.
#' @export
#'
#' @examples sampleDatacas <- EnrichData$sampleDatacas
#' enrichDB <- EnrichData$enrichDB
#' EnrichParam <- KEGGEnrichPlotPanel(sampleDatacas, enrichDB = enrichDB, pvalcutoff = 0.05, IDtype = 2, cateIdx = 1)
#' EnrichResultList <- Enrichment(EnrichParam)
#' dbChoice <- enrichDB
#' EnrichPlot(dbChoice,EnrichResultList)

EnrichPlot <- function(db, EnrichResultList){
  Table <- EnrichResultList$Table.Result
  Table <- Table[order(as.numeric(Table$P.Val),-as.numeric(Table$Count)), ]


  ## define functions
  bar.f <- NULL
  pie.f <- NULL
  if (is.null(db)==TRUE){
    return (NULL)
  }else if (db=="kegg"){
    return (NULL)
  }else if (db=="smpdb"){
    Table <- data.frame(Category.name=as.character(Table$SMPBD.name),
                        Count=as.numeric(Table$Count),
                        P.Val=as.numeric(Table$P.Val))

    bar.f <- function(tab){
      graph <- plot_bar(tab)
      print("smpbd plot out")
      return (graph)
    }
    pie.f<- function(tab){
      lbls1 <- as.character(tab$Category.name)
      lbls2 <- as.character(tab$Count)
      lbls <- paste(lbls1,lbls2,sep="\n")
      x <- as.numeric(tab$Count)
      color <- terrain.colors(length(x))
      pie(x ,labels=lbls1, col=color)
      legend("topright",lbls2 , cex=0.8, fill=color)
    }

  }else if (db=="cfam"){

    Table <- data.frame(Category.name=as.character(Table$Category.name),
                        Count=as.numeric(Table$Count),
                        P.Val=as.numeric(Table$P.Val))

    bar.f <- function(tab){
      graph <- plot_bar(tab)
      print("smpbd plot out")
      return (graph)
    }
    pie.f <- function(tab){
      lbls1 <- as.character(tab$Category.name)
      lbls2 <- as.character(tab$Count)
      lbls <- paste(lbls1,lbls2,sep="\n")
      x <- as.numeric(tab$Count)
      color <- terrain.colors(length(x))
      pie(x ,labels=lbls1, col=color,clockwise = TRUE)
      legend("topright" ,lbls2 , cex=0.8, fill=color)
    }
  }else if (db=="foodb"){

    Table <- data.frame(Category.name=as.character(Table$Category.name),
                        Count=as.numeric(Table$Count),
                        P.Val=as.numeric(Table$P.Val))

    bar.f <- function(tab){
      graph <- plot_bar(tab)
      print("smpbd plot out")
      return (graph)
    }

    pie.f <- function(tab){
      lbls1 <- as.character(tab$Category.name)
      lbls2 <- as.character(tab$Count)
      lbls <- paste(lbls1,lbls2,sep="\n")
      x <- as.numeric(tab$Count)
      color <- terrain.colors(length(x))
      pie(x ,labels=lbls1, col=color,clockwise = TRUE)
      legend("topright" ,lbls2 , cex=0.8, fill=color)
    }
  }else if (db=="biofunc"){

    Table <- data.frame(Category.name=as.character(Table$Biofunction.name),
                        Count=as.numeric(Table$Count),
                        P.Val=as.numeric(Table$P.Val))

    bar.f <- function(tab){
      graph <- plot_bar(tab)
      print("smpbd plot out")
      return (graph)
    }
    pie.f <- function(tab){
      lbls1 <- as.character(tab$Biofunction.name)
      lbls2 <- as.character(tab$Count)
      lbls <- paste(lbls1,lbls2,sep="\n")
      x <- as.numeric(tab$Count)
      color <- terrain.colors(length(x))
      pie(x ,labels=lbls1, col=color)
      legend("topright",lbls2 , cex=0.8, fill=color)
    }

  }else if (db=="tcm"){

    Table <- data.frame(Category.name=as.character(Table$Category.name),
                        Count=as.numeric(Table$Count),
                        P.Val=as.numeric(Table$P.Val))

    bar.f <- function(tab){
      graph <- plot_bar(tab)
      print("smpbd plot out")
      return (graph)
    }
    pie.f <- function(tab){
      lbls1 <- as.character(tab$Category.name)
      lbls2 <- as.character(tab$Count)
      lbls <- paste(lbls1,lbls2,sep="\n")
      x <- as.numeric(tab$Count)
      color <- terrain.colors(length(x))
      pie(x ,labels=lbls1, col=color)
      legend("topright",lbls2 , cex=0.8, fill=color)
    }
  }else if (db=="spectax"){

    Table <- data.frame(Category.name=as.character(Table$Category.name),
                        Count=as.numeric(Table$Count),
                        P.Val=as.numeric(Table$P.Val))

    bar.f <- function(tab){
      graph <- plot_bar(tab)
      print("smpbd plot out")
      return (graph)
    }
    pie.f <- function(tab){
      lbls1 <- as.character(tab$Category.name)
      lbls2 <- as.character(tab$Count)
      lbls <- paste(lbls1,lbls2,sep="\n")
      x <- as.numeric(tab$Count)
      color <- terrain.colors(length(x))
      pie(x ,labels=lbls1, col=color,clockwise = TRUE)
      legend("topright" ,lbls2 , cex=0.8, fill=color)
    }
  }else if (db=="toxin"){

    Table <- data.frame(Category.name=as.character(Table$Category.name),
                        Count=as.numeric(Table$Count),
                        P.Val=as.numeric(Table$P.Val))

    bar.f <- function(tab){
      graph <- plot_bar(tab)
      print("smpbd plot out")
      return (graph)
    }
    pie.f <- function(tab){
      lbls1 <- as.character(tab$Category.name)
      lbls2 <- as.character(tab$Count)
      lbls <- paste(lbls1,lbls2,sep="\n")
      x <- as.numeric(tab$Count)
      color <- terrain.colors(length(x))
      pie(x ,labels=lbls1, col=color,clockwise = TRUE)
      legend("topright" ,lbls2 , cex=0.8, fill=color)
    }
  }

  rowcnt <- nrow(Table)
  if (rowcnt<=3){
    return (pie.f(Table))
  }else if (rowcnt<=15){
    return (bar.f(Table))
  }else{
    return (bar.f(Table[1:15,]))
  }
}



EnrichData <- function() {
  utils::data(list="EnrichData", package="LargeMetabo")
  get("EnrichData", envir = .GlobalEnv)
}

EnrichData <- EnrichData()






