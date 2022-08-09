#' Classify New Embedded Cells 
#'
#' Classify new embedded cells using GF-ICF transformation and K-nn algorithm.
#' Existing cells are used as training set.
#' 
#' @param data list; GFICF object
#' @param classes chareachters; Classes of already existing cells in the order of they are in colnames(data$gficf).
#' @param k integer; Number of K-nn to use for classification. Odd number less than 30 are preferred.
#' @param seed integer; Initial seed to use.
#' @param method chareachters; Which space in which apply KNN (subspace or embedded). Default is embedded (i.e., umap/tsne space) because usually gives better results then the subspace PCA/NMF.
#' @return A dataframe containing cell id and predicted classes.
#' @importFrom class knn
#' 
#' @export
classify.cells = function(data,classes,k=7,seed=18051982,method="embedded")
{
  method = base::match.arg(arg = method,choices = c("subspace","embedded"),several.ok = F)
  if (is.null(data$embedded.predicted)) {stop("Please embed first new cells!")}
  set.seed(seed)
  if(!is.factor(classes)) {classes = factor(as.character(classes))}
  if (method%in%"subspace")
  {
    res = class::knn(data$pca$cells,data$pca$pred,classes,k = k,prob = T) 
  } else {
    res = class::knn(data$embedded[,c(1,2)],data$embedded.predicted[,c(1,2)],classes,k = k,prob = T)
  }
  data$embedded.predicted$predicted.class <- as.character(res)
  data$embedded.predicted$class.prob <- attr(res,"prob")
  return(data)
}

#' Embed new cells in an existing space 
#'
#' This function embed new cells in an already existing space. For now it supports only UMAP and t-UMAP. Briefly new cells are first normalized with GF-ICF method but using as ICF weigth estimated on the existing cells and than projected in the existing PCA/LSA space before to be embedded in the already existing UMAP space via umap_transform function. 
#' 
#' @param data list; GFICF object
#' @param x Matrix; UMI counts matrix of cells to embedd.
#' @param nt integer; Number of thread to use (default 2).
#' @param seed integer; Initial seed to use.
#' @param ... Additional arguments to pass to Rtsne or umap_transform call.
#' @param verbose boolean; Icrease verbosity.
#' @return The updated gficf object.
#' @import Matrix
#' @import uwot
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom Rtsne Rtsne
#' @importFrom RcppML predict.nmf
#' 
#' @export
embedNewCells = function(data,x,nt=2,seed=18051982, verbose=TRUE, ...)
{
  if(data$reduction=="tsne") {stop("Not supported with t-SNE reduction!!")}
  
  tsmessage("Gene filtering..",verbose = verbose)
  g = union(rownames(filter_genes_cell2loc_style(data = x,data$param$cell_count_cutoff,data$param$cell_percentage_cutoff2,data$param$nonz_mean_cutoff)),rownames(data$gficf))
  x = x[rownames(x)%in%g,]
  rm(g)
  
  tsmessage("Normalize counts..",verbose = verbose)
  if (data$param$normalized){
    x <- Matrix::Matrix(edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts=x),normalized.lib.sizes = T)),sparse = T)
  }
  
  x = tf(x,verbose=verbose)
  x = idf(x,w = data$w,verbose=verbose)
  x = t(l.norm(t(x),norm = "l2",verbose=verbose))
  
  if(!is.null(data$pca$odgenes) & data$pca$rescale)
  {
    x = Matrix::Matrix(data = x,sparse = T)
    x@x <- x@x*rep(data$pca$odgenes[colnames(x),'gsf'],diff(x@p))
  }
  
  if (data$pca$type=="NMF") {
    cells = colnames(x)
    x = t(RcppML::predict.nmf(w = data$pca$genes,data = x))
    rownames(x) = cells
  } else {
    x = t(x[rownames(data$pca$genes),]) %*% data$pca$genes
  }

  if(data$reduction%in%c("tumap","umap")) {
    data$embedded.predicted = as.data.frame(uwot::umap_transform(as.matrix(x),data$uwot,verbose = verbose))
    rownames(data$embedded.predicted) = rownames(x)
    colnames(data$embedded.predicted) = c("X","Y")
  }
  
  data$pca$pred = x;rm(x);gc()
  return(data)
}


