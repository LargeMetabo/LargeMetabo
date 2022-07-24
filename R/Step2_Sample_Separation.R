
#' Title Sample_Separation
#'
#' @param finalData The matrix of dataset for sample separation.
#' @param finalLabel The label of dataset for sample separation.
#' @param clusters The number of cluster for sample separation.
#' @param method The method for sample separation. The method could be HCA, KMC, PCA and SOM.
#'
#' @return The plot of sample separation.
#' @export
#'
#' @examples finalData <- MarkerData$finalData
#' finalLabel <- MarkerData$finalLabel
#' Sample_Separation(finalData, finalLabel, clusters=2, method = "HCA")


Sample_Separation <- function(finalData, finalLabel, clusters=2, method = "HCA") {


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
  library(MASS)
  library(ropls)
  library(siggenes)
  library(e1071)


  if (method == "HCA"){

    topgenes=10
    distance="euclidean"
    clust="complete"
    scale_sel="column"
    palette="Blues"
    cluster_row=clusters
    cluster_col=2
    cluster=TRUE

    sds <- rowSds(as.matrix(finalData))
    top_n <- length(sds) * topgenes / 100
    if (length(sds) > top_n) heat.data <- finalData[rev(order(sds))[1:top_n], ]
    else heat.data <- finalData

    heat.data <- t(heat.data)

    dist_method <- function(x) dist(x, method = distance)
    hclust_method <- function(x) hclust(x, method = clust)

    heatmap.plot <- d3heatmap(heat.data, scale = scale_sel, colors = palette,
                              k_row = cluster_row, k_col = cluster_col,
                              distfun = dist_method, hclustfun = hclust_method,
                              dendrogram = if (cluster) "both" else "none")
    return(heatmap.plot)


  }else if(method == "KMC"){

    iteration=10
    algorithm="Hartigan-Wong"

    data_kmeans <- finalData
    obj_kmeans <- kmeans(data_kmeans, centers = clusters, iter.max = iteration,
                         nstart = 1, algorithm = algorithm)
    km_plot <- autoplot(obj_kmeans, data = data_kmeans, label = FALSE, label.size = 5, frame = TRUE, frame.type="norm") + theme(panel.background = element_blank())

    return(km_plot)


  }else if(method == "PCA"){

    x <- as.data.frame(t(finalData))
    pcaLabel <- finalLabel
    x$label <- as.factor(pcaLabel)
    res.pca <- prcomp(x[, -ncol(x)], scale = TRUE)
    pca4_plot.out <- fviz_pca_ind(res.pca, geom = "point",
                                  habillage=x$label, addEllipses=TRUE,
                                  ellipse.level= 0.95)+ theme_minimal()
    return(pca4_plot.out)


  }else if(method == "SOM"){

    decomp <- function(x){
      fact1 <- as.integer(sqrt(x))
      while (x%%fact1!=0)
        fact1 <- fact1-1
      return (unlist(list(fact1,x%/%fact1)))
    }

    SOM_test <- function(marty.exam, marty.label, nc = 16, nct = 2) {
      library(SOMbrero)
      set.seed(1)
      my.som <- trainSOM(x.data = marty.exam, dimension = decomp(nc), nb.save = 10,
                         maxit = 2000, scaling = "none", radius.type = "letremy")
      my.sc <- superClass(my.som, k = nct)
      result <- list(som = my.som, sc = my.sc)
      return(result)
    }


    som_cluster <- clusters

    x <- t(finalData)
    rownames(x) <- 1:nrow(x)
    res.som <- SOM_test(x, Label2Digit(finalLabel), nc = som_cluster, nct = 2)
    som4_plot.out <- plot(res.som$sc, type = "lines", show.names = TRUE)

    return(som4_plot.out)


  }

}





MarkerData <- function() {
  utils::data(list="MarkerData", package="LargeMetabo")
  get("MarkerData", envir = .GlobalEnv)
}

MarkerData <- MarkerData()






