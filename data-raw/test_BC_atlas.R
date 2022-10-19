library(Matrix)

M.raw = readRDS(file = "~/work/current/BRCA_AIRC_paper/paper_git/RData/RAW.filtered.BRCA.UMI.counts.5K.umi.rds")
load("~/work/package/gficf/data/small_BC_atlas.rda")
M.raw = M.raw[,!colnames(M.raw)%in%colnames(small_BC_atlas)]

set.seed(0)
sample = sapply(strsplit(x = colnames(M.raw),split = "_",fixed = T), function(x) x[1])
u = unique(sample)
cells=NULL
for (i in 1:length(u)) {
  cells = c(cells,sample(x=colnames(M.raw)[sample%in%u[i]],size=30))
}

test_BC_atlas = M.raw[,cells]
test_BC_atlas = test_BC_atlas[rowSums(test_BC_atlas)>0,]
usethis::use_data(test_BC_atlas,overwrite = T)
