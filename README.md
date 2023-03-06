# gficf v2 - Single-cell gene set enrichment analysis

Details on ***scGSEA and scMAP*** implemented into the version 2 of GFICF package can be found in the NAR genomics and bioinformatics manuscript [Franchini et al. 2023](https://doi.org/10.1093/nargab/lqad024)  

GFICF v1 manuscript can be found here at [Gambardella et al. 2019](https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract).  

The package also includes [Phenograph](https://biorxiv.org/cgi/content/short/2022.10.24.513476v1)
[Louvain method](https://sites.google.com/site/findcommunities/)
clustering using [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy) library
from [uwot](https://github.com/jlmelville/uwot) and a naive but fast parallel implementation
of Jaccard Coefficient estimation using [RcppParallel](https://cran.r-project.org/package=RcppParallel).
The package also include data reduction with either Principal Component Analisys (PCA) or
non-negative matrix factorization [RcppML](https://github.com/zdebruine/RcppML) before to apply t-SNE or UMAP for single cell data visualization.   

**Examples & Functionality**:  
* [Install GFICF package](https://htmlpreview.github.io/?https://github.com/gambalab/gficf/blob/master/inst/doc/installation.html)  
* [Getting Started](https://htmlpreview.github.io/?https://github.com/gambalab/gficf/blob/master/inst/doc/index.html)  
* [Single-cell Gene Set Enrichement Analysis (scGSEA)](https://htmlpreview.github.io/?https://github.com/gambalab/gficf/blob/master/inst/doc/scGSEA.html)  
* [Single-cell Mapper (scMAP)](https://htmlpreview.github.io/?https://github.com/gambalab/gficf/blob/master/inst/doc/scMAP.html)  
* Batch Effect Correction (Cooming soon) 
