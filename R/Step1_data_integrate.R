
### ============================================================================
### ============================ Data Integration  =============================
### ========================== Batch Effect Removal ============================
### ============================================================================


# ~~~~~~~~~Multi file mode: begin calculation after button clicked~~~~~~~~~~

### ============================================================================
# Step 1: Integration
### ============================================================================


#' Title Integrate_Data
#'
#' @param MutileGroup Data list including multiple datasets from multiple experiments and are integrated as a comprehensive dataset.
#' @param RTTolerance1 Set the tolerance of retention time in the first step of data integration.
#' @param mzTolerance1 Set the tolerance of mass-to-charge ratio in the first step of data integration.
#' @param RTTolerance2 Set the tolerance of retention time in the second step of data integration.
#' @param mzTolerance2 Set the tolerance of mass-to-charge ratio in the second step of integration.
#'
#' @return A matrix of an integrated dataset form multiple analytical experiments is returned.
#' @export
#'
#' @examples Integrate_Data(MutileGroup, RTTolerance1 = 10, mzTolerance1 = 0.1, RTTolerance2 = 10, mzTolerance2 = 0.1)

Integrate_Data <- function(MutileGroup, RTTolerance1 = 10, mzTolerance1 = 0.1, RTTolerance2 = 10, mzTolerance2 = 0.1){


  ############################# read  mutile group ###############################

  num_Group <- length(MutileGroup)

  names(MutileGroup) <- paste("Group",1:num_Group,sep="")

  ################################ alignment #####################################

  mutile_data <- MutileGroup

  row_num <- NULL
  for (k in 1:num_Group){
    row_num <- c(row_num, nrow(mutile_data[[k]]))

  }


  mutile_align <- list()

  #for(num in 1:max(row_num)){
  for(num in 1:max(row_num)){


    max_all <- NULL
    for (n in 1:num_Group){
      max_n <- max(mutile_data[[n]][,-((ncol(mutile_data[[n]])-1):ncol(mutile_data[[n]]))])
      max_all <- c(max_all, max_n)

    }


    max_data <- mutile_data[[which.max(max_all)]]

    for(i in 1:ncol(max_data)){
      j <- match(max_all[which.max(max_all)], max_data[,i])
      if(!is.na(j)){
        loc_max <- j
      }
    }


    tentative_ref <- max_data[loc_max,]


    ############################# alignment step 1 ###############################

    all_condition <- list()

    for(m in setdiff(1:length(max_all), which.max(max_all))){

      undata_m <- mutile_data[[m]]
      condition_m <- which(abs(as.numeric(undata_m[,"mz"])-as.numeric(tentative_ref[, "mz"])) < (0.5*mzTolerance1) & abs(as.numeric(undata_m[,"rt"])-as.numeric(tentative_ref[, "rt"])) < RTTolerance1, arr.ind=TRUE)

      if(length(condition_m)!=0){
        unalign_m <- undata_m[condition_m[1],]
      }else{
        unalign_m <- 0

      }

      all_condition[[m]] <- unalign_m

    }

    all_condition[[which.max(max_all)]] <- tentative_ref



    ############################# alignment step 2 ###############################

    if(length(all_condition)!=0){

      undata_all <- all_condition


      rt_all <- list()

      for(p in 1:length(undata_all)){

        if(undata_all[[p]][1]!=0){

          rt_all[[p]] <- undata_all[[p]][,"rt"]

        }else{

          rt_all[[p]] <- NULL
        }

      }

      if(length(unlist(rt_all))%%2 == 0){

        rt_median <- unlist(rt_all)[(length(unlist(rt_all))/2+1)]

      }else{

        rt_median <- median(unlist(rt_all))

      }




      for(p in 1:length(rt_all)){

        if(length(rt_all[[p]])!=0){

          if(rt_all[[p]]==rt_median){

            ap <- p

          }

        }

      }



      unalign_s2_all <- list()

      for(p in 1:length(undata_all)){

        undata_s2_q <- undata_all[[ap]][which(rt_all[[ap]]==rt_median),]

        if (undata_all[[p]][1]!=0){

          condition_s2_q <- which(abs(as.numeric(undata_s2_q[,"mz"])-as.numeric(undata_all[[p]][, "mz"])) < (0.5*mzTolerance2) & abs(as.numeric(undata_s2_q[,"rt"])-as.numeric(undata_all[[p]][, "rt"])) < RTTolerance2, arr.ind=TRUE)

          unalign_s2_all[[p]] <- undata_all[[p]][which.min(abs(undata_all[[p]][condition_s2_q,"mz"]-as.numeric(undata_s2_q[,"mz"]))),]

        }

      }




      align <- NULL
      align <- cbind(align, paste0("Align_",num))

      for(e in 1:num_Group){

        if (e <= length(unalign_s2_all)){

          if(length(unalign_s2_all[[e]])!=0){

            new_matrix <- unalign_s2_all[[e]]
            colnames(new_matrix) <- paste0("ds", e, ".", colnames(new_matrix))

            align <- cbind(align, new_matrix )

          }else{

            new_matrix <- matrix(0, 1, ncol(mutile_data[[e]]))
            colnames(new_matrix) <- paste0("ds", e, ".", colnames(mutile_data[[e]]))

            align <- cbind(align, new_matrix)
            align <- as.data.frame(align)


          }

        }else{

          new_matrix <- matrix(0, 1, ncol(mutile_data[[e]]))
          colnames(new_matrix) <- paste0("ds", e, ".", colnames(mutile_data[[e]]))

          align <- cbind(align, new_matrix)
          align <- as.data.frame(align)


        }

      }


      colnames(align)[1] <- "Align_ID"

      mutile_align <- rbind(mutile_align, align)


      for(s in 1:length(undata_all)){


        if(s <= length(unalign_s2_all)){

          row_del <- match(row.names(unalign_s2_all[[s]]),row.names(mutile_data[[s]]))

          if(length(row_del)!=0){

            mutile_data[[s]] <- mutile_data[[s]][-which(row.names(unalign_s2_all[[s]])==row.names(mutile_data[[s]])),]
          }

        }

      }




    }else{


      align <- NULL
      align <- cbind(align, paste0("Align_",num))
      for(g in 1:num_Group){

        un_match <- match(row.names(tentative_ref),row.names(mutile_data[[g]]))

        if(length(which(!is.na(un_match)))==0){

          new_matrix <- matrix(0, 1, ncol(mutile_data[[g]]))
          colnames(new_matrix) <- paste0("ds", g, ".", colnames(mutile_data[[g]]))

          align <- cbind(align, new_matrix)
          align <- as.data.frame(align)


        }else{

          new_matrix <- tentative_ref
          colnames(new_matrix) <- paste0("ds", g, ".", colnames(new_matrix))

          align <- cbind(align, new_matrix )
          align <- as.data.frame(align)

        }

      }


      colnames(align)[1] <- "Align_ID"

      colnames(mutile_align)[1] <- "Align_ID"

      mutile_align <- rbind(mutile_align, align)


      tentative_ref

      for(h in 1:num_Group){

        row_del <- match(row.names(tentative_ref),row.names(mutile_data[[h]]))

        if(!is.na(row_del)){

          mutile_data[[h]] <- mutile_data[[h]][-which(row.names(tentative_ref)==row.names(mutile_data[[h]])),]
        }

      }

    }

  }

  mutile_align <- as.data.frame(mutile_align)

  mutile_align_all <- as.data.frame(lapply(mutile_align, as.numeric))


  row.names(mutile_align_all) <- mutile_align[, 1]
  mutile_align_all <- mutile_align_all[, -1]

  #write.csv(mutile_align_all, file="mutile_align_all.csv")
  #mutile_align_all <- read.csv("mutile_align_all.csv",header=TRUE, row.names = 1)



  return(mutile_align_all)

}


### ============================================================================
# Step 2: Batch Effect Removal
### ============================================================================

#' Title Removal_Batch
#'
#' @param MutileAlign The integrated dataset from multiple analytical experiments.
#' @param n The number of multiple analytical experiments.
#' @param algorithm The algorithm for removing batch effects in multiple analytical experiments. The algorithm can be BMC/PAMR, ComBat/EB, GlobalNorm or None.
#'
#' @return A matrix of a dataset after removing batch effects for the integrated dataset.
#' @export
#'
#' @examples Removal_Batch(MutileAlign, n = 3, algorithm = "BMC/PAMR")

Removal_Batch <- function(MutileAlign, n = 3, algorithm = "BMC/PAMR") {

  ##**************************************************************************##


  # ============================================================================ #
  # The several algorithms for batch effect remove.
  # ============================================================================ #

  #-------------------------------------------------------------------------------
  # [0] apply translog2 function for datasets.
  #-------------------------------------------------------------------------------

  translog2 <- function(x) {
    a <- lapply(x, log2)
    return(a)
  }

  #-------------------------------------------------------------------------------
  # [1] Basic function, mergNONE.
  # Notes: Combine esets without any additional transformation. Similar to 'combine' function.
  #-------------------------------------------------------------------------------

  mergeNONE <- function(x) {
    y <- x[[1]]
    for (i in 2:length(x)) {
      y <- cbind(y, x[[i]])
    }
    return(y)
  }

  #-------------------------------------------------------------------------------
  # [2] apply geneBMCESet on list of eSets
  # Notes: In they successfully applied a technique similar to z-score normalization for merging
  # breast cancer datasets. They transformed the data by batch mean-centering, which means that
  # the mean is subtracted.
  #-------------------------------------------------------------------------------

  geneBMCESet <- function(eset)
  {
    # SM: Changed / in - since exprs already on log scale
    # This corresponds to shifting mean to 0
    eset <- eset - rowMeans(eset);
    return(eset);
  }

  mergeBMC <- function(esets)
  {
    esets <- lapply(esets, geneBMCESet);
    eset <- mergeNONE(esets);
    return(eset)
  }

  #-------------------------------------------------------------------------------
  # [3] apply geneNormESet on list of eSets
  # One of the simplest mathematical transformations to make datasets more comparable
  # is z-score normalization. In this method, for each gene expression value in each study
  # separately all values are altered by subtracting the mean of the gene in that dataset
  # divided by its standard deviation
  #-------------------------------------------------------------------------------

  geneNormESet <- function(eset)
  {
    eset <- apply(eset, 1, function(x){x-mean(x, na.rm=TRUE)});
    eset <- apply(eset, 1, function(x){x/sd(x, na.rm=TRUE)});
    return(eset);
  }

  mergeGENENORM <- function(esets)
  {
    esets <- lapply(esets, geneNormESet);
    eset <-  mergeNONE(esets);
  }


  #-------------------------------------------------------------------------------
  # [4] apply ComBat method on lists of eSets.
  #-------------------------------------------------------------------------------

  bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
  postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
  postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}

  mergeCOMBAT <- function(esets){
    print("combat ininininini")
    raw_merged = esets
    batchInfo = NULL;
    for(i in 1:length(esets)) {
      #向量,有多少个平台就有多少个i，batchInfo的长度与样本的数目一样长
      batchInfo = c(batchInfo, rep(i,ncol(esets[[i]])));
    }
    sampleAll=NULL
    for (i in 1:length(esets)){
      sampleAll<- c(sampleAll,colnames(esets[[i]]))
    }

    saminfo = cbind(sampleAll,sampleAll,batchInfo)
    colnames(saminfo) = c("Array name", "Sample name", "Batch")
    dat=NULL

    for(i in 1:length(raw_merged)){
      dat=cbind(dat,as.matrix(raw_merged[[i]]))
    }
    design <- design.mat(saminfo)
    batches <- list.batch(saminfo)
    n.batch <- length(batches)
    n.batches <- sapply(batches, length)
    n.array <- sum(n.batches)
    B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))
    grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
    var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
    stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
    if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}
    s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))
    batch.design <- design[,1:n.batch]
    gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
    delta.hat <- NULL
    for (i in batches)
    {
      delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
    }
    print("half0.5")
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, aprior)
    b.prior <- apply(delta.hat, 1, bprior)

    ##Find EB batch adjustments

    gamma.star <- delta.star <- NULL
    for (i in 1:n.batch)
    {
      temp <- it.sol(s.data[,batches[[i]]], gamma.hat[i,], delta.hat[i,], gamma.bar[i], t2[i], a.prior[i], b.prior[i])
      gamma.star <- rbind(gamma.star,temp[1,])
      print("tree")
      delta.star <- rbind(delta.star,temp[2,])

    }

    bayesdata <- s.data
    j <- 1
    for (i in batches)
    {
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
      j <- j+1
    }

    bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean

    eset=raw_merged
    eset=bayesdata
    print("combat outoutouotuotuto")
    return(eset)
  }

  design.mat <- function(saminfo) {
    #第三列为batch
    tmp <- which(colnames(saminfo) == 'Batch')
    #转成因子变量
    tmp1 <- as.factor(saminfo[,tmp])
    #nlevels(tmp1)统计有多少个平台
    #msg("  => Found",nlevels(tmp1),'batches');
    design <- build.design(tmp1,start=1)
    #统计出来array,sample,batch外还有哪些注释信息
    ncov <- ncol(as.matrix(saminfo[,-c(1:2,tmp)]))
    #msg("  => Found",ncov,'covariate(s)');
    #如果还有其余的注释信息
    if(ncov>0){
      for (j in 1:ncov){
        tmp1 <- as.factor(as.matrix(saminfo[,-c(1:2,tmp)])[,j])
        design <- build.design(tmp1,des=design)
      }
    }
    design
  }

  build.design <- function(vec, des=NULL, start=2)
  {
    tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
    for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
    cbind(des,tmp)
  }

  list.batch <- function(saminfo) {
    tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'Batch')])
    batches <- NULL
    #提出每一个样本对应的平台信息，即属于哪一个批次
    for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
    batches
  }

  aprior <- function(gamma.hat) {
    m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2
  }

  it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001)
  {
    n <- apply(!is.na(sdat),1,sum)
    g.old <- g.hat
    d.old <- d.hat
    change <- 1
    count <- 0
    while(change>conv){
      g.new <- postmean(g.hat,g.bar,n,d.old,t2)
      sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
      d.new <- postvar(sum2,n,a,b)
      change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
      g.old <- g.new
      d.old <- d.new
      count <- count+1
    }
    adjust <- rbind(g.new, d.new)
    rownames(adjust) <- c("g.star","d.star")
    adjust
  }

  #-------------------------------------------------------------------------------
  # Finally, betch effect removing function was implemented.
  #-------------------------------------------------------------------------------


  BER <- function(x, method = "BMC/PAMR") {

    if (method == "BMC/PAMR")
      Merged.data <- mergeBMC(x)
    else if (method == "ComBat/EB")
      Merged.data <- mergeCOMBAT(x)
    else if (method == "GlobalNorm")
      Merged.data <- mergeGENENORM(x)
    else if (method == "Log2"){
      x <- lapply(x, log2)
      Merged.data <- mergeNONE(x)
    }
    else  if(method == "None")
      Merged.data <- mergeNONE(x)

    return(Merged.data)
  }


  # In the end.
  #-------------------------------------------------------------------------------


  ##**************************************************************************##







  esets <- list()
  p.mat <- matrix(0, n, 2)

  for (i in 1:n) {
    x <- grep(paste0("ds", i), names(MutileAlign))
    p.mat[i, 1] <- x[1]
    p.mat[i, 2] <- x[length(x)-2]
    esets[[i]] <- MutileAlign[, p.mat[i, 1]:p.mat[i, 2]]
  }

  names(esets) <- paste0("eset", 1:n)

  Merged <- BER(esets, method=algorithm)
  rownames(Merged) <- paste("Align", as.character(1:(nrow(Merged))), sep="_")

  #colnames(Merged)[1] <- "Align_ID"

  Merged[sapply(Merged, simplify = 'matrix', is.infinite)] <- 0
  Merged[sapply(Merged, simplify = 'matrix', is.na)] <- 0

  return(Merged)

}



MutileGroup <- function() {
  utils::data(list="MutileGroup", package="LargeMetabo")
  get("MutileGroup", envir = .GlobalEnv)
}

MutileGroup <- MutileGroup()






MutileAlign <- function() {
  utils::data(list="MutileAlign", package="LargeMetabo")
  get("MutileAlign", envir = .GlobalEnv)
}

MutileAlign <- MutileAlign()







