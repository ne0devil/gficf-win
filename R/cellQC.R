#' Cell QC 
#'
#' Filter Cells with low gene ratio detection and high MT ratio.
#' Loess and GAM regression are used to fit relationships between the number of UMI and either the ratio of detected genes or the MT ratio.
#' 
#' @param counts Matrix; Raw counts matrix
#' @param organism chareachters; Organism (supported Homo Sapiens and Mus musculus).
#' @param plot boolean; If regression plots must be showed.
#' @param verbose boolean; Icrease verbosity.
#' @return The updated gficf object.
#' @import AnnotationHub
#' @import ensembldb
#' 
#' @export
filter.Cells = function(counts,organism="Homo sapiens",plot=F,verbose=T) {
  organism = base::match.arg(arg = organism,choices = c("Homo sapiens","Mus musculus"),several.ok = F)
  
  metadata = data.frame(cell.id = colnames(counts),nUMI=Matrix::colSums(counts),nGenes=Matrix::colSums(counts!=0),stringsAsFactors = F)
  rownames(metadata) = metadata$cell.id
  metadata$geneRatio = metadata$nGenes/metadata$nUMI
  
  tsmessage("... Retrieving gene annotation",verbose = verbose)
  ah <- AnnotationHub()
  # Access the Ensembl database for organism
  ahDb <- query(ah,pattern = c(organism,"EnsDb"), ignore.case = TRUE)
  # Acquire the latest annotation files
  id <- tail(rownames(mcols(ahDb)),n=1)
  # Download the appropriate Ensembldb database
  edb <- ah[[id]]
  # Extract gene-level information from database
  annotations <- subset(genes(edb,return.type = "data.frame"),seq_name%in%c(as.character(1:22),"X","Y","MT"))
  # Extract IDs for mitochondrial genes
  mt = annotations$gene_id[annotations$seq_name%in%"MT"]
  # Number of UMIs assigned to mitochondrial genes
  metadata$mtUMI <- Matrix::colSums(counts[which(rownames(counts) %in% mt),], na.rm = T)
  # Calculate of mitoRatio per cell
  metadata$mitoRatio <- metadata$mtUMI/metadata$nUMI
  
  tsmessage("... Filtering Cells by Gene/nUMI ~ log(UMI)",verbose = verbose)
  tmp = plot.UMIxGene(metadata = metadata,method = "less",fdr.th = .1,plot = plot,family = "loess")
  metadata$covFilter = !tmp$toremove
  a = sum(metadata$covFilter)
  b = nrow(metadata)
  tsmessage(paste0("Cells passing the coverage filter ",a," out of ",b," (",round(a/b*100,2),")"),verbose = verbose)
  
  tsmessage("... Filtering Cells by mtRatio ~ log(UMI)")
  tmp = plot.UMIxMT(metadata,method="greater",plot = plot,family = "loess")
  metadata$mtFilter = !tmp$toremove
  a = sum(metadata$mtFilter)
  tsmessage(paste0("Cells passing the MT filter ",a," out of ",b," (",round(a/b*100,2),")"),verbose = verbose)
  
  data = list()
  data$counts = counts[,metadata$cell.id[metadata$covFilter & metadata$mtFilter]];
  rm(counts);gc()
  data$QC.metadata = metadata[colnames(data$counts),]
  data$ann.hub.id = id
  return(data)
}

#'
#' @import ggplot2
#' @importFrom MASS fitdistr
#' @importFrom mgcv gam
#' 
plot.UMIxMT = function(metadata,method="greater",fdr.th=0.1,plot=F,family="gam")
{
  if (family == "loess") {fit <- stats::loess(formula = mitoRatio ~ log(nUMI),data = metadata,span = 1, degree = 1,family = "gaussian")}
  if (family == "poly") {fit <- stats::glm(formula = mitoRatio ~ poly(log(nUMI),degree = 2,raw = T),data = metadata)}
  if (family == "gam") {fit <- mgcv::gam(mitoRatio ~ log(nUMI), data = metadata)}
  metadata$ypred = predict(fit,metadata)
  metadata$diff = metadata$mitoRatio - metadata$ypred
  fitDist = MASS::fitdistr(metadata$diff, "normal")
  
  if (method == "two.sided") {
    metadata$p.val[metadata$diff>0] = pnorm(metadata$diff[metadata$diff>0], mean = fitDist$estimate[1], sd = fitDist$estimate[2], lower.tail = FALSE)
    metadata$p.val[metadata$diff<0] = 1 - pnorm(metadata$diff[metadata$diff<0], mean = fitDist$estimate[1], sd = fitDist$estimate[2], lower.tail = FALSE)
  }
  
  if (method == "greater") {metadata$p.val = stats::pnorm(metadata$diff, mean = fitDist$estimate[1], sd = fitDist$estimate[2], lower.tail = FALSE)}
  if (method == "less") {metadata$p.val = 1 - stats::pnorm(metadata$diff, mean = fitDist$estimate[1], sd = fitDist$estimate[2], lower.tail = FALSE)}
  metadata$p.val = stats::p.adjust(metadata$p.val,method = "fdr")
  metadata$toremove = metadata$p.val<fdr.th;

  p1 = ggplot(data = metadata,aes(x=log(nUMI),y=mitoRatio)) + geom_point(aes(color=toremove),shape=19,size=.25) + theme_bw() + theme(legend.position = "none")
  if (family == "gam") {p1 = p1 + stat_smooth(method = gam, formula = y ~ s(x))}
  if (family == "loess") {p1 = p1 + stat_smooth(method = loess, formula = y ~ x,span = 1)}
  if (family == "poly") {p1 = p1 + stat_smooth(method = glm, formula = y ~ poly(x, 2, raw = TRUE))}
  if(plot) {print(p1)}
  
  return(metadata)
}

#'
#' @import ggplot2
#' @importFrom MASS fitdistr
#' @importFrom mgcv gam
#' 
plot.UMIxGene = function(metadata,method="less",fdr.th=0.1,plot=F,family="gam")
{
  #fit <- stats::glm(formula = geneRatio ~ poly(log(nUMI), d, raw = TRUE),data = metadata)
  if (family == "loess") {fit <- stats::loess(formula = geneRatio ~ log(nUMI),data = metadata,span = 1, degree = 2,family = "gaussian")}
  if (family == "gam") {fit <- mgcv::gam(geneRatio ~ log(nUMI), data = metadata)}
  metadata$ypred = predict(fit,metadata)
  metadata$diff = metadata$geneRatio - metadata$ypred
  fitDist = MASS::fitdistr(metadata$diff, "normal")
  
  if (method == "two.sided") {
    metadata$p.val[metadata$diff>0] = stats::pnorm(metadata$diff[metadata$diff>0], mean = 0, sd = fitDist$estimate[2], lower.tail = FALSE)
    metadata$p.val[metadata$diff<0] = 1 - stats::pnorm(metadata$diff[metadata$diff<0], mean = fitDist$estimate[1], sd = fitDist$estimate[2], lower.tail = FALSE)
  }
  if (method == "greater") {metadata$p.val = stats::pnorm(metadata$diff, mean = fitDist$estimate[1], sd = fitDist$estimate[2], lower.tail = FALSE)}
  if (method == "less") {metadata$p.val = 1 - stats::pnorm(metadata$diff, mean = fitDist$estimate[1], sd = fitDist$estimate[2], lower.tail = FALSE)}
  
  metadata$p.val = p.adjust(metadata$p.val,method = "fdr")
  metadata$toremove = metadata$p.val<fdr.th;
  
  p1 = ggplot(data = metadata,aes(x=log(nUMI),y=geneRatio)) + geom_point(aes(color=toremove),shape=19,size=.25) + theme_bw() + theme(legend.position = "none")
  if (family == "gam") {p1 = p1 + stat_smooth(method = gam, formula = y ~ s(x))}
  if (family == "loess") {p1 = p1 + stat_smooth(method = loess, formula = y ~ x,span = 1)}
  if(plot) {print(p1)}

  return(metadata)
}

