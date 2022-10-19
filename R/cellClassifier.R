#' Classify New Embedded Cells 
#'
#' Classify new embedded cells using GF-ICF transformation and K-nn algorithm.
#' Existing cells are used as training set.
#' 
#' @param data list; GFICF object
#' @param classes chareachters; Classes of already existing cells in the order of they are in colnames(data$gficf).
#' @param k integer; Number of K-nn to use for classification. Odd number less than 30 are preferred.
#' @param seed integer; Initial seed to use.
#' @param knn_method string; a string specifying the method. Valid methods are 'euclidean', 'manhattan', 'chebyshev', 'canberra', 'braycurtis', 'pearson_correlation', 'simple_matching_coefficient', 'minkowski' (by default the order 'p' of the minkowski parameter equals k), 'hamming', 'mahalanobis', 'jaccard_coefficient', 'Rao_coefficient'.
#' @param knn_weights_fun value; there are various ways of specifying the kernel function. NULL value (default) correspond to unweighted KNN algorithm. See the details section of KernelKnn package for more values.
#' @param nt numeric; Number of thread to use (default is 0, i.e. all possible CPUs - 1)
#' @return A dataframe containing cell id and predicted classes.
#' @importFrom KernelKnn KernelKnn
#' 
#' @export
classify.cells = function(data,classes,k=7,seed=18051982,knn_method="euclidean",knn_weights_fun=NULL,nt=0)
{
  if (is.null(data$embedded.predicted)) {stop("Please embed first new cells!")}
  if (nt==0) {nt = ifelse(detectCores()>1,detectCores()-1,1)}
  set.seed(seed)
  if(!is.factor(classes)) {classes = factor(as.character(classes))}
  
  res = KernelKnn(data = data$embedded[,c(1,2)], TEST_data = data$embedded.predicted[,c(1,2)], y = as.numeric(classes), k = k ,h = 1,method = knn_method, knn_weights_fun, threads = nt, regression = F,Levels = 1:length(levels(classes)))
  colnames(res) = levels(classes)
  res <- apply(res, 1, function(x) {i = which.max(x); return(data.frame(predicted.class=names(x)[i],class.prob=x[i],stringsAsFactors = F))} )
  res <- do.call("rbind",res)
  rownames(res) <- NULL
  #res = class::knn(data$embedded[,c(1,2)],data$embedded.predicted[,c(1,2)],classes,k = k,prob = T)
  data$embedded.predicted$predicted.class <- res$predicted.class
  data$embedded.predicted$class.prob <- res$class.prob
  return(data)
}

#' Embed new cells in an existing space 
#'
#' This function embed new cells in an already existing space. For now it supports only UMAP and t-UMAP. Briefly new cells are first normalized with GF-ICF method but using as ICF weigth estimated on the existing cells and than projected in the existing PCA/NMF space before to be embedded in the already existing UMAP space via umap_transform function. 
#' 
#' @param data list; GFICF object
#' @param x Matrix; UMI counts matrix of cells to embedd.
#' @param nt integer; Number of thread to use (default 2).
#' @param seed integer; Initial seed to use.
#' @param normalize boolean; If counts must be normalized before to be rescaled with GFICF.
#' @param verbose boolean; Icrease verbosity.
#' @return The updated gficf object.
#' @import Matrix
#' @import uwot
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom Rtsne Rtsne
#' @import RcppML
#' 
#' @export
scMAP = function(data,x,nt=2,seed=18051982, normalize=TRUE,verbose=TRUE)
{
  if(data$reduction=="tsne") {stop("Not supported with t-SNE reduction!!")}
  if(length(intersect(rownames(data$gficf),rownames(x)))==0) {stop("No common genes between the two dataset! Please check if gene identifiers beween two dataset correspond")}
  
  tsmessage("Gene filtering..",verbose = verbose)
  g = union(rownames(filter_genes_cell2loc_style(data = x,data$param$cell_count_cutoff,data$param$cell_percentage_cutoff2,data$param$nonz_mean_cutoff)),rownames(data$gficf))
  x = x[rownames(x)%in%g,]
  rm(g)
  
  if (normalize){
    tsmessage("Normalize counts..",verbose = verbose)
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


