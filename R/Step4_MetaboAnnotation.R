

#' Title MetaboAnnotation
#'
#' @param AnnotaData The value of mass-to-charge ratio for metabolite annotation.
#' @param masstole The tolerance of m/z for metabolite annotation.
#' @param toleUnit The unit of m/z for metabolite annotation. The unit could be 1 for Da, or 2 for ppm.
#' @param annotaDB The database of m/z for metabolite annotation. The database could be metlin or hmdb.
#' @param ionMode The mode of m/z for metabolite annotation. The mode could be pos, neg or neu.
#'
#' @return The list of results of metabolite annotation for different charges.
#' @export
#'
#' @examples AnnotaMS <- AnnotaData$AnnotaMS
#' MetaboA_res <- MetaboAnnotation(AnnotaMS)
#' MetaboA_res$`M+H-2H2O`

MetaboAnnotation <- function(AnnotaData, masstole = 10, toleUnit = 1, annotaDB = "metlin", ionMode  = "pos"){
  
  param <- list()
  
  # inputId="masstole",label="Tolerance",value=10,step=0.0001,min=0
  param$masstole <- masstole
  
  if (toleUnit==1){
    param$toleUnit <- "Da"
  }else if (toleUnit==2){
    param$toleUnit <- "ppm"
  }else{
  }
  
  # inputId="annotaDB", label="Annotation MS_database:", choices=c("METLIN"="metlin", "HMDB"="hmdb")
  
  param$annotaDB <- annotaDB
  
  # inputId="ionMode",  label="Mode:", choices=c("Positive"="pos", "Negative"="neg", "Neutral"="neu")
  
  param$ionMode  <- ionMode
  
  
  #param$charge <- charge
  
  mzlist <- as.numeric(AnnotaData)
  param$mzlist <- mzlist[!is.na(mzlist)]
  
  AnnotaParam <- param
  
  
  ####
  MPchargeName <- c("M+H-2H2O","M+H-H2O","M-H","M-H2O+NH4","M+H","M+Li",
                    "M+NH4","M+Na","M+CH3OH+H","M+K","M+ACN+H","M+2Na-H",
                    "M+ACN+Na","M+2H","M+H+Na","M+2Na","M+3H","M+2H+Na","M+H+2Na")
  MNchargeName <- c("M-H2O-H","M-H","M+F","M+Na-2H","M+Cl","M+K-2H","M+FA-H","M+CH3COO","M-2H","M-3H")
  
  HPchargeName <- c("M+3H" ,"M+2H+Na"  ,"M+H+2Na"  ,"M+2H","M+H+NH4","M+H+Na","M+H+K"    ,"M+ACN+2H",
                    "M+2Na","M+2ACN+2H","M+3ACN+2H","M+H" ,"M+NH4"  ,"M+Na"  ,"M+CH3OH+H","M+ACN+H" ,
                    "M+IsoProp+H","M+ACN+Na","M+DMSO+H","M+2Na-H","M+IsoProp+Na+H","2M+H","2M+NH4",
                    "2M+Na","2M+3ACN+2H","2M+K","2M+3H2O+H","2M+ACN+Na")
  
  HNchargeName <- c("M-2H","M-3H","M-H2O-H","M-H","M+Na-2H","M+Cl","M+K-2H","M+FA-H",
                    "M+Hac-H","M+Br","M+TFA-H","2M-H","2M+FA-H","2M+Hac-H","3M-H")
  
  
  if (annotaDB=="metlin" && ionMode=="pos"){
    charge <- 1:19
    param$charge <- charge
    AnnotaParam$charge <- charge
    chargeName <- MPchargeName[charge]
  }else if(annotaDB=="metlin" && ionMode=="neg"){
    charge <- 1:10
    param$charge <- charge
    AnnotaParam$charge <- charge
    chargeName <- MNchargeName[AnnotaParam$charge]
  }else if(annotaDB=="hmdb" && ionMode=="pos"){
    charge <- 1:19
    param$charge <- charge
    AnnotaParam$charge <- charge
    chargeName <- HPchargeName[AnnotaParam$charge]
  }else if(annotaDB=="hmdb" && ionMode=="neg"){
    charge <- 1:15
    param$charge <- charge
    AnnotaParam$charge <- charge
    chargeName <- HNchargeName[AnnotaParam$charge]
  }else if (ionMode=="neu"){
    chargeName <-  "M"
  }else{
    chargeName <-  NULL
  }
  
  
  ####
  #load("D:/Code/Rpackage/LargeMetabo/data/annotation_data.RData")
  
  annotation_data <- function() {
    utils::data(list="annotation_data", package="LargeMetabo")
    get("annotation_data", envir = .GlobalEnv)
  }
  
  annotation_data <- annotation_data()
  
  
  
  Annotation <- function(param){
    
    adduct_list <- annotation_data$adduct_list
    
    ionList <- adduct_list[intersect(grep(param$annotaDB, adduct_list$compound_db),grep(param$ionMode, adduct_list$ion_mode)),]
    
    calcMZLowerBound <- NULL
    calcMZUpperBound <- NULL
    if (param$ionMode == "pos" || param$ionMode == "neg"){
      calcMZLowerBound <- function(mz, unit, db, tole, charge, aum, mcof){
        eps <- tole * ifelse(unit=="ppm", mz*1e-6, 1)
        return (charge* (mz-eps-aum) / ifelse(db=="hmdb",mcof,1) )
      }
      calcMZUpperBound <- function(mz, unit, db, tole, charge, aum, mcof){
        eps <- tole * ifelse(unit=="ppm", mz*1e-6, 1)
        return (charge* (mz+eps-aum) / ifelse(db=="hmdb",mcof,1) )
      }
    }else if (param$ionMode == "neu"){
      calcMZLowerBound <- function(mz, unit, db, tole, charge, aum, mcof){
        if (unit=="ppm")
          return (mz/(1+tole*1e-6))
        else
          return (mz-tole)
      }
      calcMZUpperBound <- function(mz, unit, db, tole, charge, aum, mcof){
        if (unit=="ppm")
          return (mz/(1-tole*1e-6))
        else
          return (mz+tole)
      }
    }
    
    unit <- param$toleUnit
    db   <- param$annotaDB
    tole <- param$masstole
    
    AnnotationTables <- list()
    
    for (addIdx in 1:length(param$charge)){
      adduct <- param$charge[addIdx]
      queryResult <- NULL
      for (mz in param$mzlist){
        minMz <- calcMZLowerBound(mz,unit,db, tole, ionList$charge[adduct], ionList$aum[adduct], ionList$mcof[adduct])
        maxMz <- calcMZUpperBound(mz,unit,db, tole, ionList$charge[adduct], ionList$aum[adduct], ionList$mcof[adduct])
        if(db=='metlin'){
          
          MS_database <- annotation_data$MS_database
          queryLines <- MS_database[which(MS_database$MetlinMass>=minMz & MS_database$MetlinMass<=maxMz),c("LargeMetabo_ID", "MetlinMass", "MID", "Common_name", "Annotation")]
          
        }else{
          queryLines <- MS_database[which(MS_database$MonoisotopicMass>=minMz & MS_database$MonoisotopicMass<=maxMz & MS_database$Source==param$annotaDB),c("LargeMetabo_ID", "MonoisotopicMass", "MID", "Common_name", "Annotation")]
          
        }
        
        queryLinesWithQMZ <-  cbind(QueryMZ=rep(mz,nrow(queryLines)), queryLines)
        queryResult <- rbind(queryResult, queryLinesWithQMZ)
        queryResult <- unique(queryResult)
      }
      AnnotationTables[addIdx] <- list(queryResult)
    }
    
    if (param$ionMode=="neu"){
      names(AnnotationTables) <- c("M")
    }else{
      names(AnnotationTables) <- ionList$name[param$charge]
    }
    
    return (AnnotationTables)
  }
  
  
  ##---------------------------------------------------------------------------------------------------------
  
  parameter <- AnnotaParam
  resultTable <- NULL
  annotationResult <- Annotation(parameter)
  
  return(annotationResult)
  
  
}










#' Title AnnotaTandom
#'
#' @param Parent_mass The value of parent ion mass for metabolite annotation of tandom mass spectrum.
#' @param TandomData The value of MS/MS peak list (m/z & Intensity) for metabolite annotation of tandom mass spectrum.
#' @param massTandom The tolerance of parent ion mass for metabolite annotation of tandom mass spectrum.
#' @param toleUnitTandom The unit of parent ion mass for metabolite annotation of tandom mass spectrum. The unit could be 1 for Da, or 2 for ppm.
#' @param massmzTandom The tolerance of MS/MS peak mass for metabolite annotation of tandom mass spectrum.
#' @param toleUnitmzTandom The unit of MS/MS peak mass for metabolite annotation of tandom mass spectrum. The unit could be 1 for Da, or 2 for ppm.
#' @param ModeTandom The ionization mode for metabolite annotation of tandom mass spectrum. The mode could be Positive or Negative.
#' @param ionEnergy The CID energy for metabolite annotation of tandom mass spectrum. The mode could be low(10V) Medium(25V) High(40V) or All.
#'
#' @return The parameters of metabolite annotation of tandom mass spectrum.
#' @export
#'
#' @examples Parent_mass <- AnnotaData$Parent_mass
#' TandomData <- AnnotaData$TandomData
#' AnnotaTandom(Parent_mass, TandomData)

AnnotaTandom <- function(Parent_mass = 181.04, TandomData, massTandom = 0.1, toleUnitTandom = 1, massmzTandom = 0.5, toleUnitmzTandom = 1, ModeTandom = "Positive", ionEnergy = "low(10V)"){
  
  ####
  paramTandom <- list()
  
  # Parent_mass,label="Parent Ion Mass (Da)",value=181.04,step=0.0001,min=0
  
  paramTandom$Parent_mass <- Parent_mass
  
  # massTandom,label="Tolerance of Parent Ion",value=0.1,step=0.0001,min=0
  
  paramTandom$massTandom <- massTandom
  
  # toleUnitTandom,label="Unit",choices=list("Da"=1,"ppm"=2),selected=1,inline=T
  
  if (toleUnitTandom==1){
    paramTandom$toleUnitTandom <- "Da"
  }else if (toleUnitTandom==2){
    paramTandom$toleUnitTandom <- "ppm"
  }
  
  # massmzTandom",label="Tolerance of Mass/Charge",value=0.5,step=0.0001,min=0
  
  paramTandom$massmzTandom <- massmzTandom
  
  # "toleUnitmzTandom",label="Unit",choices=list("Da"=1,"ppm"=2),selected=1,inline=T
  
  if (toleUnitmzTandom==1){
    paramTandom$toleUnitmzTandom <- "Da"
  }else if (toleUnitmzTandom==2){
    paramTandom$toleUnitmzTandom <- "ppm"
  }
  
  #
  #param$annotaDB <- input$annotaDB
  #param$ionMode  <- input$ionMode
  # ModeTandom, label="Ionization Mode", choices=c("Positive"="Positive", "Negative"="Negative")))
  
  if (ModeTandom=="Positive"){
    paramTandom$ModeTandom <- "Positive"
  }else if(ModeTandom=="Negative"){
    paramTandom$ModeTandom <- "Negative"
  }else{
    paramTandom$ModeTandom <- NULL
  }
  
  # inputId="ionEnergy",  label="CID Energy",  choices=c("low(10V)"="low(10V)", "Medium(25V)"="Medium(25V)", "High(40V)"="High(40V)", "All"="All")))
  
  if (ionEnergy=="low(10V)"){
    paramTandom$ionEnergy <- "low(10V)"
  }else if(ionEnergy=="Medium(25V)"){
    paramTandom$ionEnergy <- "Medium(25V)"
  }else if(ionEnergy=="High(40V)"){
    paramTandom$ionEnergy <- "High(40V)"
  }else if(ionEnergy=="All"){
    paramTandom$ionEnergy <- "All"
  }else{
    paramTandom$ionEnergy <- NULL
  }
  
  # Include_predicted", label=strong("Include predicted spectra"), value = FALSE, width = NULL
  
  paramTandom$Include_predicted <- FALSE
  
  # TandomData
  colnames(TandomData) <- c("V1","RI")
  
  mzTandomlist <- as.data.frame(TandomData)
  paramTandom$mzTandomlist <- mzTandomlist[!is.na(mzTandomlist)]
  
  AnnotaParamTandom <- paramTandom
  
  return(AnnotaParamTandom)
  
}

AnnotaTandem <- function(Parent_mass = 181.04, TandomData, massTandom = 0.1, toleUnitTandom = 1, massmzTandom = 0.5, toleUnitmzTandom = 1, ModeTandom = "Positive", ionEnergy = "low(10V)"){
  
  ####
  paramTandom <- list()
  
  # Parent_mass,label="Parent Ion Mass (Da)",value=181.04,step=0.0001,min=0
  
  paramTandom$Parent_mass <- Parent_mass
  
  # massTandom,label="Tolerance of Parent Ion",value=0.1,step=0.0001,min=0
  
  paramTandom$massTandom <- massTandom
  
  # toleUnitTandom,label="Unit",choices=list("Da"=1,"ppm"=2),selected=1,inline=T
  
  if (toleUnitTandom==1){
    paramTandom$toleUnitTandom <- "Da"
  }else if (toleUnitTandom==2){
    paramTandom$toleUnitTandom <- "ppm"
  }
  
  # massmzTandom",label="Tolerance of Mass/Charge",value=0.5,step=0.0001,min=0
  
  paramTandom$massmzTandom <- massmzTandom
  
  # "toleUnitmzTandom",label="Unit",choices=list("Da"=1,"ppm"=2),selected=1,inline=T
  
  if (toleUnitmzTandom==1){
    paramTandom$toleUnitmzTandom <- "Da"
  }else if (toleUnitmzTandom==2){
    paramTandom$toleUnitmzTandom <- "ppm"
  }
  
  #
  #param$annotaDB <- input$annotaDB
  #param$ionMode  <- input$ionMode
  # ModeTandom, label="Ionization Mode", choices=c("Positive"="Positive", "Negative"="Negative")))
  
  if (ModeTandom=="Positive"){
    paramTandom$ModeTandom <- "Positive"
  }else if(ModeTandom=="Negative"){
    paramTandom$ModeTandom <- "Negative"
  }else{
    paramTandom$ModeTandom <- NULL
  }
  
  # inputId="ionEnergy",  label="CID Energy",  choices=c("low(10V)"="low(10V)", "Medium(25V)"="Medium(25V)", "High(40V)"="High(40V)", "All"="All")))
  
  if (ionEnergy=="low(10V)"){
    paramTandom$ionEnergy <- "low(10V)"
  }else if(ionEnergy=="Medium(25V)"){
    paramTandom$ionEnergy <- "Medium(25V)"
  }else if(ionEnergy=="High(40V)"){
    paramTandom$ionEnergy <- "High(40V)"
  }else if(ionEnergy=="All"){
    paramTandom$ionEnergy <- "All"
  }else{
    paramTandom$ionEnergy <- NULL
  }
  
  # Include_predicted", label=strong("Include predicted spectra"), value = FALSE, width = NULL
  
  paramTandom$Include_predicted <- FALSE
  
  # TandomData
  colnames(TandomData) <- c("V1","RI")
  
  mzTandomlist <- as.data.frame(TandomData)
  paramTandom$mzTandomlist <- mzTandomlist[!is.na(mzTandomlist)]
  
  AnnotaParamTandom <- paramTandom
  
  return(AnnotaParamTandom)
  
}




##########
AnnotationTandom <- function(AnnotaParamTandom){
  
  calcMZLowerTandom <- NULL
  calcMZLowerTandom <- function(mz, unit, tole){
    if (unit=="ppm")
      return (mz/(1+tole*1e-6))
    else
      return (mz-tole)
  }
  
  calcMZUpperTandom <- NULL
  calcMZUpperTandom <- function(mz, unit, tole){
    if (unit=="ppm")
      return (mz/(1-tole*1e-6))
    else
      return (mz+tole)
  }
  
  mzT <- AnnotaParamTandom$Parent_mass
  unitT <- AnnotaParamTandom$toleUnitTandom
  toleT <- AnnotaParamTandom$massTandom
  
  #
  queryResult <- NULL
  minMzT <- calcMZLowerTandom(mzT, unitT, toleT)
  maxMzT <- calcMZUpperTandom(mzT, unitT, toleT)
  
  #
  #load("D:/Code/Rpackage/LargeMetabo/data/annotationTandom_data.RData")
  
  annotationTandom_data01 <- function() {
    utils::data(list="annotationTandom_data01", package="LargeMetabo")
    get("annotationTandom_data01", envir = .GlobalEnv)
  }
  
  annotationTandom_data01 <- annotationTandom_data01()
  
  
  annotationTandom_data02 <- function() {
    utils::data(list="annotationTandom_data02", package="LargeMetabo")
    get("annotationTandom_data02", envir = .GlobalEnv)
  }
  
  annotationTandom_data02 <- annotationTandom_data02()
  
  
  annotationTandom_data <- rbind(annotationTandom_data01, annotationTandom_data02)
  
  
  Annotation_db <- annotationTandom_data
  options(digits=2)
  AnnotationTandom <- Annotation_db[which(Annotation_db$Mass>=minMzT & Annotation_db$Mass<=maxMzT & Annotation_db$Energy == AnnotaParamTandom$ionEnergy & Annotation_db$Ionization == AnnotaParamTandom$ModeTandom),c("LargeMetabo_ID", "Label", "Spectrum", "Mass", "Name", "Annotation")]
  
  return(AnnotationTandom)
  
}


##########
#' Title annotaDataTandom
#'
#' @param AnnotaParamTandom The parameters of metabolite annotation of tandom mass spectrum.
#'
#' @return The table of results for metabolite annotation of tandom mass spectrum.
#' @export
#'
#' @examples Parent_mass <- AnnotaData$Parent_mass
#' TandomData <- AnnotaData$TandomData
#' AnnotaParamTandom <- AnnotaTandom(Parent_mass, TandomData)
#' annotaDataTandom(AnnotaParamTandom)

annotaDataTandom <- function(AnnotaParamTandom){
  
  ####
  AnnotationTables <- AnnotationTandom(AnnotaParamTandom)
  
  ##
  res_MSMS <- NULL
  for(i in 1:dim(AnnotationTables)[1]){
    
    res_co <- NULL
    spectrumTandom <- AnnotationTables$Spectrum[i]
    
    s1 <- unlist(strsplit(spectrumTandom,";"))
    
    s1 <- gsub(":",",",s1)
    
    res_s1 <- NULL
    for(i in 1:length(s1)){
      
      str_a <- NULL
      str_s1 <- strsplit(s1[i], ",")
      str_a <- unlist(str_s1)
      res_s1 <- c(res_s1, str_a)
      
      
    }
    
    s11 <- as.numeric(unlist(res_s1))
    
    #
    s2 <- matrix(s11,length(s1),2,byrow = TRUE)
    
    colnames(s2) <- c("V1","RI")
    
    ###
    mergeTolerance <- function(x, y, tolerance = 1e-5) {
      colnames(x) <-
        c("V1", 2:ncol(x)) #suppresses error warning 'duplicate column names'
      colnames(y) <-
        c("V1", 2:ncol(y)) #suppresses error warning 'duplicate column names'
      mrg <- merge(x, y, by = "V1", all = TRUE)
      mrg[is.na(mrg)] <- 0
      i <- 1
      while (!is.na(mrg[(i + 1), 1])) {
        if (abs(mrg[i, 1] - mrg[(i + 1), 1]) <= mrg[i, 1] * tolerance) {
          mrg[i, 1] <- (mrg[i, 1] + mrg[(i + 1), 1]) / 2
          mrg[i,-1] <- mrg[i,-1] + mrg[(i + 1),-1]
          mrg <- mrg[-(i + 1),]
          i <- i + 1
          colnames(mrg) <-
            c("V1", 2:ncol(mrg))
          ##suppresses error warning 'duplicate column names'
        } else {
          i <- i + 1
        }
      }
      mrg
    }
    
    cossim <- function(x, y, type = c("spectrum", "neutral_losses"),
                       mzTolerance = 1e-5) {
      colnames(x) <- NULL
      colnames(y) <- NULL
      mm <- mergeTolerance(x, y, tolerance = mzTolerance)
      sum(sqrt(mm[, 2]) * sqrt(mm[, 3])) /
        (sqrt(sum(mm[, 2])) * sqrt(sum(mm[, 3])))
    }
    ###
    
    
    res_co <- cossim(TandomData, s2)
    
    res_co <- round(res_co,2)
    
    res_MSMS <- rbind(res_MSMS, res_co)
    
  }
  
  colnames(res_MSMS) <- "Fit"
  row.names(res_MSMS) <- 1:dim(res_MSMS)[1]
  
  AnnotationTandom <- rbind(NULL,NULL)
  for (i in 1:dim(AnnotationTables)[1]){
    
    bi_a <- cbind(AnnotationTables[i,],res_MSMS[i])
    
    AnnotationTandom <- rbind(AnnotationTandom, bi_a)
    
  }
  
  colnames(AnnotationTandom) <- c(colnames(AnnotationTables),"Fit")
  
  AnnotationTandom1 <- AnnotationTandom[,-which(colnames(AnnotationTandom)=="Spectrum")]
  
  Include_predicted <- AnnotaParamTandom$Include_predicted
  
  
  if (Include_predicted==FALSE){
    
    AnnotationTandom2 <- AnnotationTandom1[,-which(colnames(AnnotationTandom1)=="Label")]
    
  }else{
    AnnotationTandom2 <- AnnotationTandom1[AnnotationTandom1[,"Label"]=="experimental",-which(colnames(AnnotationTandom1)=="Label")]
  }
  
  annotationResultTandom <- AnnotationTandom2
  
  ####
  mass <- AnnotaParamTandom$mzTandomlist
  
  table_s <- annotationResultTandom
  
  annotaDataTandom <- table_s[rev(order(table_s$Fit)),]
  
  return(annotaDataTandom)
  
}



annotaDataTandem <- function(AnnotaParamTandom){
  
  ####
  AnnotationTables <- AnnotationTandom(AnnotaParamTandom)
  
  ##
  res_MSMS <- NULL
  for(i in 1:dim(AnnotationTables)[1]){
    
    res_co <- NULL
    spectrumTandom <- AnnotationTables$Spectrum[i]
    
    s1 <- unlist(strsplit(spectrumTandom,";"))
    
    s1 <- gsub(":",",",s1)
    
    res_s1 <- NULL
    for(i in 1:length(s1)){
      
      str_a <- NULL
      str_s1 <- strsplit(s1[i], ",")
      str_a <- unlist(str_s1)
      res_s1 <- c(res_s1, str_a)
      
      
    }
    
    s11 <- as.numeric(unlist(res_s1))
    
    #
    s2 <- matrix(s11,length(s1),2,byrow = TRUE)
    
    colnames(s2) <- c("V1","RI")
    
    ###
    mergeTolerance <- function(x, y, tolerance = 1e-5) {
      colnames(x) <-
        c("V1", 2:ncol(x)) #suppresses error warning 'duplicate column names'
      colnames(y) <-
        c("V1", 2:ncol(y)) #suppresses error warning 'duplicate column names'
      mrg <- merge(x, y, by = "V1", all = TRUE)
      mrg[is.na(mrg)] <- 0
      i <- 1
      while (!is.na(mrg[(i + 1), 1])) {
        if (abs(mrg[i, 1] - mrg[(i + 1), 1]) <= mrg[i, 1] * tolerance) {
          mrg[i, 1] <- (mrg[i, 1] + mrg[(i + 1), 1]) / 2
          mrg[i,-1] <- mrg[i,-1] + mrg[(i + 1),-1]
          mrg <- mrg[-(i + 1),]
          i <- i + 1
          colnames(mrg) <-
            c("V1", 2:ncol(mrg))
          ##suppresses error warning 'duplicate column names'
        } else {
          i <- i + 1
        }
      }
      mrg
    }
    
    cossim <- function(x, y, type = c("spectrum", "neutral_losses"),
                       mzTolerance = 1e-5) {
      colnames(x) <- NULL
      colnames(y) <- NULL
      mm <- mergeTolerance(x, y, tolerance = mzTolerance)
      sum(sqrt(mm[, 2]) * sqrt(mm[, 3])) /
        (sqrt(sum(mm[, 2])) * sqrt(sum(mm[, 3])))
    }
    ###
    
    
    res_co <- cossim(TandomData, s2)
    
    res_co <- round(res_co,2)
    
    res_MSMS <- rbind(res_MSMS, res_co)
    
  }
  
  colnames(res_MSMS) <- "Fit"
  row.names(res_MSMS) <- 1:dim(res_MSMS)[1]
  
  AnnotationTandom <- rbind(NULL,NULL)
  for (i in 1:dim(AnnotationTables)[1]){
    
    bi_a <- cbind(AnnotationTables[i,],res_MSMS[i])
    
    AnnotationTandom <- rbind(AnnotationTandom, bi_a)
    
  }
  
  colnames(AnnotationTandom) <- c(colnames(AnnotationTables),"Fit")
  
  AnnotationTandom1 <- AnnotationTandom[,-which(colnames(AnnotationTandom)=="Spectrum")]
  
  Include_predicted <- AnnotaParamTandom$Include_predicted
  
  
  if (Include_predicted==FALSE){
    
    AnnotationTandom2 <- AnnotationTandom1[,-which(colnames(AnnotationTandom1)=="Label")]
    
  }else{
    AnnotationTandom2 <- AnnotationTandom1[AnnotationTandom1[,"Label"]=="experimental",-which(colnames(AnnotationTandom1)=="Label")]
  }
  
  annotationResultTandom <- AnnotationTandom2
  
  ####
  mass <- AnnotaParamTandom$mzTandomlist
  
  table_s <- annotationResultTandom
  
  annotaDataTandom <- table_s[rev(order(table_s$Fit)),]
  
  return(annotaDataTandom)
  
}


################################################################################

specplot2 <- function(n, o, list) {
  plot(
    x = list[[n]][, 1],
    y = list[[n]][, 2] / max(list[[n]][, 2]),
    col = "blue",
    type = "h",
    xlim = c((min(c(
      list[[n]][, 1], list[[o]][, 1]
    )) * 0.9), (max(c(
      list[[n]][, 1], list[[o]][, 1]
    )) * 1.1)),
    xaxs = "i",
    xlab = expression(italic(m / z)),
    ylim = c(-1.2, 1.2),
    yaxs = "i",
    yaxt = "n",
    ylab = "intensity relative to base peak",
    main = paste(names(list[n]), "Mirror plot of two spectra (red: the input, blue: the database)", names(list[o]))
  )
  points(x = list[[o]][, 1],
         y = -(list[[o]][, 2] / max(list[[o]][, 2])),
         col = "red",
         type = "h")
  abline(a = 0, b = 0)
  axis(2,
       at = seq(-1, 1, 0.5),
       labels = c(1.0, 0.5, 0.0, 0.5, 1.0))
  text(
    x = (list[[n]][, 1])[(list[[n]][, 2] / max(list[[n]][, 2])) > 0.1],
    y = (list[[n]][, 2] / max(list[[n]][, 2]))[(list[[n]][, 2] / max(list[[n]][, 2])) > 0.1],
    labels = round((list[[n]][, 1])[(list[[n]][, 2] / max(list[[n]][, 2])) > 0.1], 4),
    pos = 3,
    cex = 0.75
  )
  text(
    x = (list[[o]][, 1])[(list[[o]][, 2] / max(list[[o]][, 2])) > 0.1],
    y = -((list[[o]][, 2] / max(list[[o]][, 2]))[(list[[o]][, 2] / max(list[[o]][, 2])) > 0.1]),
    labels = round((list[[o]][, 1])[(list[[o]][, 2] / max(list[[o]][, 2])) > 0.1], 4),
    pos = 1,
    cex = 0.75
  )
}


##########################################

#' Title annotaTandom_plot
#'
#' @param AnnotaParamTandom The parameters of metabolite annotation of tandom mass spectrum.
#' @param TandomData The value of MS/MS peak list (m/z & Intensity) for metabolite annotation of tandom mass spectrum.
#'
#' @return The plot of results for metabolite annotation of tandom mass spectrum.
#' @export
#'
#' @examples Parent_mass <- AnnotaData$Parent_mass
#' TandomData <- AnnotaData$TandomData
#' AnnotaParamTandom <- AnnotaTandom(Parent_mass, TandomData)
#' annotaDataTandom(AnnotaParamTandom)
#' annotaTandom_plot(AnnotaParamTandom, TandomData)

annotaTandom_plot <- function(AnnotaParamTandom, TandomData){
  
  mass <- AnnotaParamTandom$mzTandomlist
  
  annotResTan <- AnnotationTandom(AnnotaParamTandom)
  
  AnnotationTables <- annotResTan
  plot_cand <- AnnotationTables$LargeMetabo_ID
  
  shortlist <- list()
  shortlist[[1]] <- TandomData
  
  AnnotationTables <- annotResTan
  
  MSMSid <- plot_cand
  
  mSpectrum <- AnnotationTables[which(AnnotationTables[,"LargeMetabo_ID"]==MSMSid),"Spectrum"]
  
  
  
  m1 <- unlist(strsplit(mSpectrum,";"))
  
  m1 <- gsub(":",",",m1)
  
  res_m1 <- NULL
  for(i in 1:length(m1)){
    
    str_am <- NULL
    str_m1 <- strsplit(m1[i], ",")
    str_am <- unlist(str_m1)
    res_m1 <- c(res_m1, str_am)
    
    
  }
  
  m11 <- as.numeric(unlist(res_m1))
  
  #
  annotatedSpeclist2 <- matrix(m11,length(m1),2,byrow = TRUE)
  
  colnames(annotatedSpeclist2) <- c("V1","RI")
  
  shortlist[[2]] <- annotatedSpeclist2
  
  specplot2(1,2,shortlist)
  
}




annotaTandem_plot <- function(AnnotaParamTandom, TandomData){
  
  mass <- AnnotaParamTandom$mzTandomlist
  
  annotResTan <- AnnotationTandom(AnnotaParamTandom)
  
  AnnotationTables <- annotResTan
  plot_cand <- AnnotationTables$LargeMetabo_ID
  
  shortlist <- list()
  shortlist[[1]] <- TandomData
  
  AnnotationTables <- annotResTan
  
  MSMSid <- plot_cand
  
  mSpectrum <- AnnotationTables[which(AnnotationTables[,"LargeMetabo_ID"]==MSMSid),"Spectrum"]
  
  
  
  m1 <- unlist(strsplit(mSpectrum,";"))
  
  m1 <- gsub(":",",",m1)
  
  res_m1 <- NULL
  for(i in 1:length(m1)){
    
    str_am <- NULL
    str_m1 <- strsplit(m1[i], ",")
    str_am <- unlist(str_m1)
    res_m1 <- c(res_m1, str_am)
    
    
  }
  
  m11 <- as.numeric(unlist(res_m1))
  
  #
  annotatedSpeclist2 <- matrix(m11,length(m1),2,byrow = TRUE)
  
  colnames(annotatedSpeclist2) <- c("V1","RI")
  
  shortlist[[2]] <- annotatedSpeclist2
  
  specplot2(1,2,shortlist)
  
}


AnnotaData <- function() {
  utils::data(list="AnnotaData", package="LargeMetabo")
  get("AnnotaData", envir = .GlobalEnv)
}

AnnotaData <- AnnotaData()






