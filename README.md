# gficf v2 package overview

An R implementation of the 
Gene Frequency - Inverse Cell Frequency method for single cell data
normalization [(Gambardella et al. 2019)](https://www.frontiersin.org/articles/10.3389/fgene.2019.00734/abstract).
The package also includes [Phenograph](https://www.cell.com/cell/fulltext/S0092-8674(15)00637-6)
[Louvain method](https://sites.google.com/site/findcommunities/)
clustering using [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy) library
from [uwot](https://github.com/jlmelville/uwot) and a naive but fast parallel implementation
of Jaccard Coefficient estimation using [RcppParallel](https://cran.r-project.org/package=RcppParallel).
The package also include data reduction with either Principal Component Analisys (PCA) or
non-negative matrix factorization [RcppML](https://github.com/zdebruine/RcppML) before to apply t-SNE or UMAP for single cell data visualization.   

**Examples & Functionality**:  
* [Getting Started](https://htmlpreview.github.io/?https://github.com/gambalab/gficf/blob/master/inst/doc/index.html)  
* [Single-cell Gene Set Enrichement Analysis (scGSEA)](https://htmlpreview.github.io/?https://github.com/gambalab/gficf/blob/master/inst/doc/scGSEA.html)  
* [Single-cell Mapper (scMAP)](https://htmlpreview.github.io/?https://github.com/gambalab/gficf/blob/master/inst/doc/scMAP.html)  
* Batch Effect Correction (Cooming soon) 


# Installation


#### 1. OS required Steps (Officially supported only Linux)

`gficf` makes use of `Rcpp`, `RcppParallel` and `RcppGSL`. So you have to carry out
a few extra steps before being able to build this package. The steps are reported below for each platform.


##### 1.1 Ubuntu/Debian

You need gsl dev library to successfully install RcppGSL library.
On Ubuntu/Debian systems this can be accomplished by runnuing the command `sudo apt-get install libgsl-dev libcurl4-openssl-dev libssl-dev libxml2-dev` from the terminal.


##### 1.2 Mac OS X

1.2.1 Open terminal and run `xcode-select --install` to install the command line developer tools.

1.2.1. We than need to install gsl libraries. This can be done via [Homebrew](https://brew.sh/). So, still from terminal
```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```
and than use `homebrew` to install gsl with following command
```bash
brew install gsl
```


##### 1.3 Windows

1.3.1 Skip this first step if you are using RStudio because it will ask you automatically. Otherwise install  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and ensure  `path\to\Rtools\bin` is on your path.   

1.3.2 [Download gsl library for Windows](https://sourceforge.net/projects/gnu-scientific-library-windows/) from sourceforge and exctract it in `C:\` or where you want.   

1.3.3 Open R/Rstudio and before to istall the package from github exec the following command in the R terminal.
```R
# Change the path if you installed gsl librarie not in the default path.
# Be sure to use the format '"path/to/gsl-xxx_mingw-xxx/gsl-xxx-static"'
# In this way " characters will be mainteined and spaces in the path preserved if there are.

# For example for gsl-2.2.1 compiled with mingw-6.2.0:
Sys.setenv(GSL_LIBS = '"C:/gsl-2.2.1_mingw-6.2.0/gsl-2.2.1-static"')
```


#### 2. After above OS specific steps

Exec in R terminal the following commands
```R
# Install required bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
BiocManager::install(setdiff(c("sva","edgeR", "fgsea"),rownames(installed.packages())),update = F)

if(!require(devtools)){ install.packages("devtools")}
devtools::install_github("gambalab/gficf")
```
