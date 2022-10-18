library(gficf)
library(ggplot2)
library(ggbeeswarm)
library(plyr)

M.raw = readRDS(file = "~/work/current/BRCA_AIRC_paper/paper_git/RData/RAW.filtered.BRCA.UMI.counts.5K.umi.rds")

set.seed(0)
sample = sapply(strsplit(x = colnames(M.raw),split = "_",fixed = T), function(x) x[1])
u = unique(sample)
cells=NULL
for (i in 1:length(u)) {
  if (u[i]=="HDQP1") {
    cells = c(cells,sample(x=colnames(M.raw)[sample%in%u[i]],size=110))
  } else {
    cells = c(cells,sample(x=colnames(M.raw)[sample%in%u[i]],size=150))
  }
}

small_BC_atlas = M.raw[,cells]
small_BC_atlas = small_BC_atlas[rowSums(small_BC_atlas)>0,]
usethis::use_data(small_BC_atlas,overwrite = T)
