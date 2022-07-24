

### Metohd-1. Fold Change -----------------------------------------------------#

FC_test <- function(data, labels, method = median) {

  if (missing(data))
    stop("** FAILURE : 'data' is missing **")
  if (missing(labels))
    stop("** FAILURE : 'labels' is missing **")

  data1.ave <- apply(data[, which(labels == 0)], 1, method, na.rm = TRUE)
  data2.ave <- apply(data[, which(labels == 1)], 1, method, na.rm = TRUE)
  fold.change <- data2.ave - data1.ave

  result <- data.frame(Variable = names(fold.change), ExpreInControl = data1.ave,
                       ExpreInCase = data2.ave, FoldChange = round(fold.change, 2))
  rownames(result) <- NULL
  return(result)


}

FC_barplot <- function(mat, lab, num, titl) {

  x <- mat[1:num, ]
  y <- ifelse(lab=="0", "control", "case")
  colnames(x) <- y
  case.pos <- colnames(x)=="case"
  control.pos <- colnames(x)=="control"
  case_data <- x[, case.pos]
  control_data <- x[, control.pos]

  data_list <- list()

  for (i in 1:num) {

    js <- as.numeric(control_data[i, ])
    os <- as.numeric(case_data[i, ])

    data_list[[i*2-1]] <- js
    data_list[[i*2]] <- os

  }

  names(data_list) <- paste(rep(rownames(x), each = 2), rep(c("Control", "Case"), num), sep = "_")
  #op <- par(mar= c(10, 3, 1, 2))
  cols <- rep(c("green", "red"), times = num)
  boxplot(data_list, las = 2, col = cols, notch =TRUE, main = titl)
  #par(op)
}


### Metohd-2. PLS-DA -----------------------------------------------------#

library(mixOmics)

PLSDA_test <- function(mat, label, ncom = 2, ite = 500, tole = 1e-06, nearZV = TRUE) {

  X <- t(mat)
  Y <- label
  plsda.res <- plsda(X, Y, ncomp = ncom, max.iter = ite, tol = tole, near.zero.var = nearZV)
  col.res <- as.numeric(as.factor(Y))
  plotIndiv(plsda.res, ind.names = TRUE, add.legend =TRUE, plot.star = TRUE, plot.centroid = TRUE,
            abline.line = TRUE, plot.ellipse = TRUE)
}


PLSDA.VIP <- function (model, graph = FALSE)
{
  if (packageVersion("mixOmics") < "5.0.2") {
    stop(paste("you must update 'mixOmics' to version >= 5.0.2 (actual: ",
               packageVersion("mixOmics"), ")", sep = ""))
  }
  VIP <- mixOmics::vip(model)
  tab <- as.data.frame(VIP[order(VIP[, ncol(VIP)], decreasing = TRUE),
                           ncol(VIP)])
  colnames(tab) <- "VIP"
  if (graph) {
    opar <- par()
    on.exit(suppressWarnings(par(opar)))
    par(mar = c(5, 8, 2, 2), las = 1)
    g <- barplot(rev(tab$VIP), horiz = TRUE,
                 xlab = paste("VIP (", ncol(VIP), ifelse(ncol(VIP) > 1, " axes)", " axis)"),
                              sep = ""))
    mtext(rev(rownames(tab)), side = 2, line = 1, at = g,
          cex = 0.7)
    abline(h = g, lty = 3, col = "grey40")
    abline(v = 1, lty = 2, col = "red")
  }
  result <- list(tab = tab, sup1 = rownames(tab)[which(tab$VIP > 1)])
  class(result) <- "PLSDA.VIP"
  return(result)
}


### Metohd-3. OPLS-DA -----------------------------------------------------#

library(ropls)
OPLSDA_test <- function(mat, label, cross = 7, per = 0, scale = "standard") {
  X <- t(mat)
  Y <- as.factor(label)
  oplsda <- opls(X, Y, predI = 1, orthoI = NA, plotL = TRUE, crossvalI = cross, permI = per, scaleC = scale)
  res <- oplsda@vipVn
  cpds <- data.frame(Variable = names(res), VIP = res)
  cpds <- cpds[rev(order(cpds$VIP)), ]
  rownames(cpds) <- NULL

  oplsda_res <- list()
  oplsda_res$oplsda <- oplsda
  oplsda_res$VIP <- cpds
  return(oplsda_res)
}

### Metohd-4. T Test -----------------------------------------------------#

FWER.Bonf <-  function (pval){
  m <- length(pval)
  signif <- rep(FALSE,m)

  ## Bonferroni
  resFDR<-mt.rawp2adjp(pval, proc=c("Bonferroni"))
  adjp <- resFDR$adjp[order(resFDR$index),]

  return (adjp[,"Bonferroni"])
}

FDR.BH<-  function (pval, TST=FALSE,q){
  m <- length(pval)
  signif <- rep(FALSE,m)

  resFDR<-mt.rawp2adjp(pval, proc=c("BH"))
  adjp <- resFDR$adjp[order(resFDR$index),]

  if (TST){

    r<-mt.reject(adjp[,"BH"],(q/(1+q)))$r
    outFDR<-adjp[,"BH"]*((m-r)/m)
  }
  else{
    outFDR<-adjp[,"BH"]
  }

  return (outFDR)
}

multiple.correction <- function(pval, typeFDR, q) {
  if (missing(pval))
    stop("Error : raw pvalues")
  if (missing(typeFDR))
    stop("Error : no typeFDR")
  if (!typeFDR %in% c("FWER", "FDR-BH", "FDR-TST", "qvalue"))
    stop("Error : unknown typeFDR")
  print(paste("typeFDR=", typeFDR))
  if (typeFDR == "FWER")
    padj <- FWER.Bonf(pval)
  else if (typeFDR == "FDR-BH")
    padj <- FDR.BH(pval)
  else if (typeFDR == "FDR-TST")
    padj <- FDR.BH(pval, TST = TRUE, q = 0.05)
  else if (typeFDR == "qvalue") {
    pi0 <- pi0.est(pval)$p0
    padj <- qvalue.cal(pval, pi0, version = 2)
  }
  return(padj)
}


T_test <- function(data, labels, typeFDR = "FDR-BH", algo = "t", q = 0.05, plot = TRUE) {
  require(multtest)
  list.exp <- colnames(data)
  list.probe <- rownames(data)
  if (is.null(list.probe)) {
    list.probe <- 1:nrow(data)
  }
  message("Launch ", algo, " test")
  test <- mt.teststat(data, labels, test = algo)
  if (algo == "t.equalvar") {
    s1 <- apply(data[, which(labels == 0)], 1, function(x) {
      return(length(which(!is.na(x))))
    })
    s2 <- apply(data[, which(labels == 1)], 1, function(x) {
      return(length(which(!is.na(x))))
    })
    ddl <- s1 + s2 - 2
  }
  else if (algo == "pairt") {
    s1 <- apply(data[, which(labels == 0)], 1, function(x) {
      return(length(which(!is.na(x))))
    })
    ddl <- s1 - 1
  }
  else if (algo == "t") {
    s1 <- apply(data[, which(labels == 0)], 1, function(x) {
      return(length(which(!is.na(x))))
    })
    s2 <- apply(data[, which(labels == 1)], 1, function(x) {
      return(length(which(!is.na(x))))
    })
    w1 <- apply(data[, which(labels == 0)], 1, var, na.rm = TRUE)/s1
    w2 <- apply(data[, which(labels == 1)], 1, var, na.rm = TRUE)/s2
    ddl <- (w1 + w2)^2/((w1^2/(s1 - 1)) + (w2^2/(s2 - 1)))
  }
  message("Calculate pval")
  pval.test <- 2 * (pt(-abs(test), ddl))
  message("Adjusted pval")
  padj <- multiple.correction(pval.test, typeFDR)
  indexsignif <- which(padj <= q)
  message(length(indexsignif), "significant genes")
  out <- data.frame(list.probe, test, pval.test, padj)
  colnames(out) <- c("probeID", "Stat", "RawpValue", "AdjpValue")
  if (plot) {
    col <- rep("black", length(test))
    col[indexsignif] <- "red"
    qqnorm(test, col = col)
    qqline(test)
  }
  return(out)
}

### Metohd-5. Chi.squared test -----------------------------------------------------#

CHIS_test <- function(x, y) {
  data <- as.data.frame(t(x))
  data$label <- y
  weights <- chi.squared(label~., data)
  chis_res <- data.frame(attr_importance = weights$attr_importance)
  rownames(chis_res) <- rownames(weights)
  return(chis_res)
}

### Metohd-6. Correction-based Method -----------------------------------------------------#

CFS_test <- function(data, label, method = "spearman") {
  cor.mat <- matrix(0, nrow(data), 2)
  for (i in 1:nrow(data)) {
    cor.mat[i, 1] <- row.names(data)[i]
    cor.mat[i, 2] <- cor(as.numeric(data[i, ]), label, method  = method)
  }
  cor.mat<-cor.mat[rev(order(cor.mat[,2])),]
  cor.mat <- as.data.frame(cor.mat)
  names(cor.mat) <- c("Variable", "Correlation coefficient")
  return(cor.mat)
}

### Metohd-7. ENTROPY_test Method -----------------------------------------------------#

ENTROPY_test <- function(mat, label) {
  x <- t(mat)
  cl <- label
  comb <- as.data.frame(cbind(x, cl))
  ga <- information.gain(cl~., comb)
  un <- symmetrical.uncertainty(cl~., comb)
  ga.order <- rev(order(ga$attr_importance))
  un.order <- rev(order(un$attr_importance))
  EBM_res <- data.frame(InfGain = rownames(ga)[ga.order], InfGain.Order = ga$attr_importance[ga.order],
                        SymmeUn = rownames(un)[un.order], SymmeUn.Order = un$attr_importance[un.order])
  return(EBM_res)
}

### Metohd-8. Linear Model and Emperical eBayes method -----------------------------------------------------#

LMEB_test <- function(data, labels) {
  if (class(data) != "matrix") data <- as.matrix(data)
  design <- cbind(Grp1 = 1, Grp2vsGrp1 = labels)
  library(limma)
  fit <- lmFit(data, design)
  fit <- eBayes(fit)
  topTable(fit, coef = 2, n = nrow(fit))
}

### Metohd-9. RELIEF test method -----------------------------------------------------#

RELIEF_test <- function(mat, label, nc = 5, ss = 20) {
  x <- t(mat)
  cl <- label
  comb <- as.data.frame(cbind(x, cl))

  weights <- relief(cl~., comb, neighbours.count = nc, sample.size = ss)
  we.order <- rev(order(weights$attr_importance))

  REL_res <- data.frame(Index = paste("Index", as.character(1:nrow(mat)), sep="_"),
                        Feature.Order = rownames(weights)[we.order],
                        Importance.Order = weights$attr_importance[we.order])
  REL_res <- as.data.frame(REL_res)
  return(REL_res)
}

### Metohd-10. RF test method -----------------------------------------------------#

RF_test <- function(mat, label, nt = 500, ntI = 300, vdf = 0.2, seed = 1) {
  set.seed(seed)
  x <- t(mat)
  cl <- factor(label)
  rf.vs1 <- varSelRF(x, cl, ntree = nt, ntreeIterat = ntI, vars.drop.frac = vdf)
  return(rf.vs1)

}

### Metohd-11. SAM test method -----------------------------------------------------#

require(siggenes)
SAM_test <- function (data, labels, nbpermut = 500, q = 0.05, plot = TRUE, method = "d.stat", var.equal = TRUE,
                      include.zero = FALSE, paired = FALSE, seed = 123) {
  if (paired == FALSE) {
    labels <- as.integer(factor(labels)) - 1
  }

  message("SAM analysis is running...")
  if (method == "d.stat"){
    output <- sam(data, cl = labels, B = nbpermut, method = "d.stat",
                  var.equal = var.equal, include.zero = include.zero,
                  rand = seed, control = samControl(delta = seq(0.1, 10, 0.05), lambda = 0.5), med = TRUE) #alter delta from 0.1 to 0.01
  }


  if (method == "wilc.stat")
    output <- sam(data, cl = labels, method = "wilc.stat")
  if (method == "chisq.stat")
    output <- sam(data, cl = labels, method = "chisq.stat", B = nbpermut, rand = seed)
  if (method == "trend.stat")
    output <- sam(data, cl = labels, method = "trend.stat", B = nbpermut, rand = seed)

  message("Create table for all genes...")
  if (length(unique(labels)) == 2 & (method == "d.stat" || method == "wilc.stat")) {
    alllist <- data.frame(rownames(data), output@d, output@p.value, output@fold)
    colnames(alllist) <- c("probeID", "Stat", "RawpValue", "FoldChange")
  }
  else {
    alllist <- data.frame(rownames(data), output@d, output@p.value)
    colnames(alllist) <- c("probeID", "Stat", "RawpValue")
  }

  alllist$probeID = as.character(alllist$probeID)
  print(output@mat.fdr[which(output@mat.fdr[, "FDR"] > 0),
                       c("Delta", "p0", "False", "Called", "FDR")])
  message("Find delta...")

  outdelta <- try(findDelta(output, fdr = q))

  if (class(outdelta) == "matrix") {
    delta <- outdelta[2, 1]
    print(paste("Delta : ", delta, sep = ""))
    message("Create SAM plot...")
    if (plot) {
      plot(output, delta)
    }
    alllist$Significant <- rep(FALSE, nrow(alllist))
    if (outdelta[2, 1] != 0) {
      delta.sum <- summary(output, delta)
      ds.mat.sig <- delta.sum@mat.sig
      rows <- ds.mat.sig$Row
      print(paste("Find", length(rows), "significant genes ..."))
      alllist[rows, "Significant"] <- TRUE
    }
  }

  message("The adjusted pvalue is a qvalue, and is not related to the significant genes found with the SAM's FDR")
  return(alllist)
}

### Metohd-12. SVM-RFE method -----------------------------------------------------#

library(e1071)
svmrfeFeatureRanking = function(x, y, cValue = 50, Scaled = FALSE){
  n = ncol(x)

  survivingFeaturesIndexes = seq(1:n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n

  while(length(survivingFeaturesIndexes)>0){
    svmModel <- svm(x[, survivingFeaturesIndexes],
                    y,
                    cost = cValue,
                    cachesize=500,
                    scale=Scaled,
                    type="C-classification",
                    kernel="linear")
    w <- t(svmModel$coefs) %*% svmModel$SV
    rankingCriteria <- w * w
    ranking <- sort(rankingCriteria, index.return = TRUE)$ix
    featureRankedList[rankedFeatureIndex] <- survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex <- rankedFeatureIndex - 1
    (survivingFeaturesIndexes <- survivingFeaturesIndexes[-ranking[1]])
  }
  return (featureRankedList)
}

### Metohd-13. WRST method -----------------------------------------------------#

WRST <- function(data, label, type = "approximate", perm = 99) {
  con.loc <- which(label == "0")
  case.loc <- which(label == "1")
  f <- function(x) {
    wil.res <- wilcox.test(x[con.loc], x[case.loc])
    wil.res$p.value
  }

  f_w <- function(x) {
    x1 <- as.vector(x[con.loc])
    y1 <- as.vector(x[case.loc])
    DV <- c(x1, y1)
    IV <- factor(rep(c("control", "case"), c(length(x1), length(y1))))

    if (type == "approximate")
      w_per <- pvalue(wilcox_test(DV ~ IV, alternative="two.sided", distribution=approximate(B=perm)))
    if (type == "exact")
      w_per <- pvalue(wilcox_test(DV ~ IV, alternative="two.sided", distribution = "exact", conf.int = TRUE))
    if (type == "asymptotic")
      w_per <- pvalue(wilcox_test(DV ~ IV, alternative="two.sided", distribution = "asymptotic", conf.int = TRUE))

    w_per
  }

  w.result <- data.frame(Variable = rownames(data), P.value = apply(data, 1, f))
  w.result
}





#################################################################################


#' Title
#'
#' @param finalData The matrix of dataset for marker identification.
#' @param finalLabel The label of dataset for marker identification.
#' @param method The method for marker identification. The method could be FC, PLS-DA, OPLS-DA, t-test, CHIS, CFS, ENTROPY, LMEB, RELIEF, RF, SAM, SVMRFE or WRST.
#'
#' @return The result of marker identification.
#' @export
#'
#' @examples finalData <- MarkerData$finalData
#' finalLabel <- MarkerData$finalLabel
#' Marker_Identify(finalData, finalLabel, method = "FC")

Marker_Identify <- function(finalData, finalLabel, method = "FC") {

  library(genefilter)
  library(d3heatmap)
  library(ggplot2)
  library(ggfortify)
  library(factoextra)
  library(igraph)
  library(FSelector)
  library(varSelRF)
  library(corrplot)
  library(mixOmics)
  library(ropls)
  library(siggenes)
  library(e1071)


  if (method == "FC"){

    marker_res <- list()

    topNum <- 10
    tailNum <- 10

    par(mfrow = c(2, 1))

    #
    FC.out <- FC_test(t(finalData), labels = finalLabel, mean)
    headFinalData <- finalData[rev(order(FC.out$FoldChange)), ]
    FC_plot.out_head <- FC_barplot(headFinalData, finalLabel, as.numeric(topNum),"boxplots of the top 10 most up-regulated markers")
    marker_res$FC_plot_head <- FC_plot.out_head

    #
    tailFinalData <- finalData[(order(FC.out$FoldChange)), ]
    FC_plot.out_tail <- FC_barplot(tailFinalData, finalLabel, as.numeric(tailNum),"boxplots of the top 10 most down-regulated markers")
    marker_res$FC_plot_tail <- FC_plot.out_tail

    par(mfrow = c(1, 1))


    #
    FC.out <- FC.out[rev(order(FC.out$FoldChange)), ]
    rownames(FC.out) <- NULL
    FC.res <- data.frame(FC.out$Variable, round(FC.out$ExpreInControl, 6), round(FC.out$ExpreInCase, 6), round(FC.out$FoldChange, 2))
    names(FC.res) <- c("Variable", "Expression Value in Control", "Expression Value in Case", "Fold Change")
    marker_res$FC_table <- FC.res

    return(marker_res)


  }else if(method == "PLS-DA"){

    prcom_num = 2
    ite_times = 500
    tol_ite = 1e-6

    #
    plsda_plot1 <- PLSDA_test(finalData, finalLabel, ncom = prcom_num, ite = ite_times, tole = as.numeric(tol_ite))

    #
    tmp <- plsda(t(finalData), finalLabel, ncomp = prcom_num, max.iter = ite_times, tol = as.numeric(tol_ite))
    res <- PLSDA.VIP(tmp, graph = FALSE)

    obj <- res$tab
    new.obj <- data.frame(Variable = rownames(obj), VIP = obj[, 1])
    rownames(new.obj) <- NULL
    VIP_value <- formatC(new.obj$VIP, format = "f", digits = 6)
    plsda.res <- data.frame(variable = new.obj$Variable, VIP = VIP_value)

    return(plsda.res)


  }else if(method == "OPLS-DA"){

    cross_num = 7
    perm_num = 0
    scaleC = "standard"

    oplsda_ret <- OPLSDA_test(finalData, finalLabel, cross = as.numeric(cross_num), per = perm_num, scale = scaleC)
    VIP_table <- oplsda_ret$VIP

    return(VIP_table)


  }else if(method == "t-test"){


    Comptype = "unequal variances"
    pvaluetype = 0.05

    #plot
    corrected_type <- c("t", "t.equalvar", "pairt")
    names(corrected_type) <- c("unequal variances", "equal variances", "paired")
    TTest.out <- T_test(finalData, finalLabel, typeFDR = "FDR-BH", algo = as.vector(corrected_type[Comptype]), plot = FALSE)
    new.TTest.out <- TTest.out[order(TTest.out$AdjpValue), ]
    colnames(new.TTest.out) <- c("Variable", "T Statistics", "P Value", "Adjusted P Value")
    rownames(new.TTest.out) <- NULL

    cols <- rep("cyan4", nrow(TTest.out))
    cols[TTest.out[, 4] < as.numeric(pvaluetype)] <- "red"
    plot(-log(TTest.out[, 4], 10), xlab = "Metabolites (AlignID start from 0)", ylab = "-log10(adjPvalue)",
         pch = 19, col = cols, main = "")
    grid()
    abline( h= -log(as.numeric(pvaluetype), 10), col = "darkorchid", lty = 2, lwd = 2)


    #table

    return(new.TTest.out)



  }else if(method == "CHIS"){

    #
    CHIS.out <- CHIS_test(finalData, finalLabel)
    cols <- rep("cyan4", nrow(CHIS.out))
    cols[CHIS.out$attr_importance > 0.2] <- "red"
    CHIS_plot <- plot(CHIS.out$attr_importance, xlab = "Metabolites (AlignID start from 0)", ylab = "Attribute Importance",
                      pch = 19, col = cols, main = "")

    #
    CHIS.res <- data.frame(Variable = rownames(CHIS.out), Importance = CHIS.out[, 1])
    CHIS.res <- CHIS.res[rev(order(CHIS.res$Importance)), ]
    rownames(CHIS.res) <- NULL

    return(CHIS.res)


  }else if(method == "CFS"){

    method="spearman"

    #
    M <- cor(finalData, method = method)
    cor_plot <- corrplot(M, order = "hclust", addrect = 2, col = rainbow(10))

    #
    CFS.out <- CFS_test(finalData, finalLabel, method = method)

    return(CFS.out)


  }else if(method == "ENTROPY"){

    res <- ENTROPY_test(finalData, finalLabel)

    return(res)


  }else if(method == "LMEB"){

    #
    LMEB.out <- LMEB_test(finalData, finalLabel)
    LMEB.out$threshold <- ifelse(abs(LMEB.out$logFC) >= 1 & LMEB.out$P.Value < 0.05, 2, 4)
    plot(LMEB.out$logFC, -log10(LMEB.out$P.Value), col = LMEB.out$threshold, pch = 20,
         ylab = expression(-log[10](P.Value)), xlab = expression(log[2](FC)))
    abline(h=2, lty = 2, col = "lightgreen")
    abline(v=c(-1, 1), lty = 2, col = "lightgreen")

    #
    return(LMEB.out)


  }else if(method == "RELIEF"){

    nc_value=5
    ss_value=20

    res <- RELIEF_test(finalData, finalLabel, nc = nc_value, ss = ss_value)

    return(res)


  }else if(method == "RF"){

    tree_num=500
    ntI_value=300
    vdf_value=0.2
    seed_value=1


    rf.vs <- RF_test(finalData, finalLabel, nt = tree_num, ntI = ntI_value, vdf = vdf_value, seed = seed_value)
    rf.rfe_res <- rf.vs$firstForest$importance

    rf.rfe_res1 <- data.frame(Variable = rownames(rf.rfe_res), control = rf.rfe_res[, 1],
                             case = rf.rfe_res[, 2], MeanDecreaseAccurary = rf.rfe_res[, 3],
                             MeanDecreaseGini = rf.rfe_res[, 4])
    rownames(rf.rfe_res1) <- NULL

    return(rf.rfe_res1)


  }else if(method == "SAM"){

    SAM.out <- SAM_test(finalData, finalLabel)

    return(SAM.out)


  }else if(method == "SVMRFE"){

    c_value = 50
    scaling_svm = FALSE

    SVMRFE.out <- svmrfeFeatureRanking(t(finalData), finalLabel, cValue = c_value, Scaled = scaling_svm)
    result_table <- data.frame(Index.of.Feature = SVMRFE.out,
                               Feature.Name = rownames(finalData)[SVMRFE.out])
    return(result_table)


  }else if(method == "WRST"){

    permType="approximate"
    permtimes=99

    WRST.out <- WRST(finalData, finalLabel, type = permType , perm = as.numeric(permtimes))
    WRST.out <- WRST.out[order(WRST.out$P.value), ]
    rownames(WRST.out) <- NULL

    return(WRST.out)


  }


}






#################################################################################



#' Title
#'
#' @param finalData The matrix of dataset for marker identification.
#' @param finalLabel The label of dataset for marker identification.
#' @param method The method for marker identification. The method could be FC, PLS-DA, OPLS-DA, t-test, CHIS, CFS, ENTROPY, LMEB, RELIEF, RF, SAM, SVMRFE or WRST.
#'
#' @return The result of marker identification.
#' @export
#'
#' @examples finalData <- MarkerData$finalData
#' finalLabel <- MarkerData$finalLabel
#' Marker_Identify(finalData, finalLabel, method = "FC")

















#' Title
#'
#' @param finalData The matrix of dataset for marker identification.
#' @param finalLabel The label of dataset for marker identification.
#' @param method The method for marker identification. The method could be FC, PLS-DA, OPLS-DA, t-test, CHIS, CFS, ENTROPY, LMEB, RELIEF, RF, SAM, SVMRFE or WRST.
#'
#' @return The ROC curve and AUC value in the SVM classification using markers.
#' @export
#'
#' @examples finalData <- MarkerData$finalData
#' finalLabel <- MarkerData$finalLabel
#' Marker_Assess(finalData, finalLabel, method = "FC")
#' 
Marker_Assess <- function(finalData, finalLabel, method = "FC") {
  
  set.seed(1)
  
  library(genefilter)
  library(d3heatmap)
  library(ggplot2)
  library(ggfortify)
  library(factoextra)
  library(igraph)
  library(FSelector)
  library(varSelRF)
  library(corrplot)
  library(mixOmics)
  library(ropls)
  library(siggenes)
  library(e1071)
  library(sampling)
  library(AUC)
  
  
  row.names(finalData) <- 1:nrow(finalData)
  
  
  if (method == "FC"){
    
    #
    FC.out <- FC_test(t(finalData), labels = finalLabel, mean)
    
    markers <- which((abs(FC.out$FoldChange) > 0.5)==TRUE)
    
    
  }else if(method == "PLS-DA"){
    
    prcom_num = 2
    ite_times = 500
    tol_ite = 1e-6
    
    #
    tmp <- plsda(t(finalData), finalLabel, ncomp = prcom_num, max.iter = ite_times, tol = as.numeric(tol_ite))
    res <- PLSDA.VIP(tmp, graph = FALSE)
    
    markers <- row.names(res$tab)[which((res$tab > 1)==TRUE)]
    markers <- as.numeric(markers)
    
    
  }else if(method == "OPLS-DA"){
    
    cross = 7
    per = 0
    scale = "standard"
    
    X <- t(finalData)
    Y <- as.factor(finalLabel)
    oplsda <- opls(X, Y, predI = 1, orthoI = NA, fig.pdfC = FALSE, crossvalI = cross, permI = per, scaleC = scale)
    res <- oplsda@vipVn
    cpds <- data.frame(Variable = names(res), VIP = res)
    cpds <- cpds[rev(order(cpds$VIP)), ]
    rownames(cpds) <- NULL
    
    VIP_table <- cpds
    
    markers <- VIP_table$Variable[which((VIP_table$VIP > 1)==TRUE)]
    markers <- as.numeric(markers)
    
    
  }else if(method == "t-test"){
    
    Comptype = "unequal variances"
    pvaluetype = 0.05
    
    #
    corrected_type <- c("t", "t.equalvar", "pairt")
    names(corrected_type) <- c("unequal variances", "equal variances", "paired")
    TTest.out <- T_test(finalData, finalLabel, typeFDR = "FDR-BH", algo = as.vector(corrected_type[Comptype]), plot = FALSE)
    new.TTest.out <- TTest.out[order(TTest.out$AdjpValue), ]
    colnames(new.TTest.out) <- c("Variable", "T Statistics", "P Value", "Adjusted P Value")
    rownames(new.TTest.out) <- NULL
    
    markers <- new.TTest.out$Variable[which((new.TTest.out$`Adjusted P Value` < 0.05)==TRUE)]
    markers <- as.numeric(markers)
    
    
  }else if(method == "CHIS"){
    
    #
    CHIS.out <- CHIS_test(finalData, finalLabel)
    
    CHIS.res <- data.frame(Variable = rownames(CHIS.out), Importance = CHIS.out[, 1])
    CHIS.res <- CHIS.res[rev(order(CHIS.res$Importance)), ]
    
    markers <- CHIS.res$Variable[which((CHIS.res$Importance > 0)==TRUE)]
    markers <- as.numeric(markers)
    
    
  }else if(method == "CFS"){
    
    method="spearman"
    
    #
    CFS.out <- CFS_test(finalData, finalLabel, method = method)
    
    markers <- as.numeric(CFS.out$Variable)[which((abs(as.numeric(CFS.out$`Correlation coefficient`)) > 0.6)==TRUE)]
    
    
  }else if(method == "ENTROPY"){
    
    res <- ENTROPY_test(finalData, finalLabel)
    
    markers <- res$InfGain[which((res$InfGain.Order > 0)==TRUE)]
    markers <- as.numeric(markers)
    
  }else if(method == "LMEB"){
    
    #
    LMEB.out <- LMEB_test(finalData, finalLabel)
    LMEB.out$threshold <- ifelse(abs(LMEB.out$logFC) >= 1 & LMEB.out$P.Value < 0.05, 2, 4)
    
    markers <- row.names(LMEB.out)[abs(LMEB.out$logFC) >= 1 & LMEB.out$adj.P.Val < 0.05]
    markers <- as.numeric(markers)
    
    
  }else if(method == "RELIEF"){
    
    nc = 5
    ss = 20
    
    x <- t(finalData)
    cl <- finalLabel
    comb <- as.data.frame(cbind(x, cl))
    
    weights <- relief(cl~., comb, neighbours.count = nc, sample.size = ss)
    
    markers <- cutoff.k(weights, 100)
    markers <- as.numeric(markers)
    
    
  }else if(method == "RF"){
    
    nt = 500
    ntI = 300
    vdf = 0.2
    seed = 1
    
    x <- t(finalData)
    cl <- factor(finalLabel)
    rf.vs1 <- varSelRF(x, cl, ntree = nt, ntreeIterat = ntI, vars.drop.frac = vdf)
    
    markers <- rf.vs1$selected.vars
    markers <- as.numeric(markers)
    #rf.vs1$firstForest$importance
    
    
  }else if(method == "SAM"){
    
    SAM.out <- SAM_test(finalData, finalLabel)
    
    markers <- SAM.out$probeID[abs(SAM.out$FoldChange) >= 1 & SAM.out$RawpValue < 0.05]
    markers <- as.numeric(markers)
    
    
  }else if(method == "SVMRFE"){
    
    c_value = 50
    scaling_svm = FALSE
    
    SVMRFE.out <- svmrfeFeatureRanking(t(finalData), finalLabel, cValue = c_value, Scaled = scaling_svm)
    result_table <- data.frame(Index.of.Feature = SVMRFE.out,
                               Feature.Name = rownames(finalData)[SVMRFE.out])
    
    markers <- result_table$Feature.Name[1:100]
    markers <- as.numeric(markers)
    
    
  }else if(method == "WRST"){
    
    permType="approximate"
    permtimes=99
    
    WRST.out <- WRST(finalData, finalLabel, type = permType , perm = as.numeric(permtimes))
    
    markers <- WRST.out$Variable[which((WRST.out$P.value < 0.05)==TRUE)]
    markers <- as.numeric(markers)
    
    
  }
  
  
  
  X_matrix <- finalData[markers, ]
  y_label <- as.factor(finalLabel)
  
  X <- t(X_matrix)
  y <- y_label
  
  # cross-validated SVM-probability plot
  
  test_data <- as.data.frame(cbind(y, X))
  colnames(test_data)[1] <- "Group"
  
  folds <- 2
  test.fold <- list()
  
  #test.fold <- split(sample(1:length(number_labels)), 1:folds)
  test.fold1 <- strata(test_data, "Group",size=(as.numeric(table(test_data$Group))/2),method="srswor")[,2]
  
  test.fold[[1]] <- test.fold1
  test.fold[[2]] <- (1:nrow(test_data))[-test.fold1]
  
  #
  all.pred.tables <-  lapply(1:folds, function(i) {
    test <- test.fold[[i]]
    Xtrain <- X[-test, ]
    ytrain <- as.factor(y[-test])
    sm <- svm(Xtrain, ytrain, prob = TRUE) # some tuning may be needed
    prob.benign <- attr(predict(sm, X[test,], prob = TRUE), "probabilities")[, 2]
    data.frame(ytest = y[test], ypred = prob.benign) # returning this
  })
  
  full.pred.table <- do.call(rbind, all.pred.tables)
  
  plot(roc(full.pred.table[, 2], full.pred.table[, 1]), col = "red")
  
  auc.value <- auc(roc(full.pred.table[, 2], full.pred.table[, 1]))
  
  auc.value <- formatC(auc.value, format = "f", digits = 3)
  
  auc_p <- paste0("The AUC of the ROC curve is ", auc.value)
  
  return(auc_p)
  
}


















MarkerData <- function() {
  utils::data(list="MarkerData", package="LargeMetabo")
  get("MarkerData", envir = .GlobalEnv)
}

MarkerData <- MarkerData()














