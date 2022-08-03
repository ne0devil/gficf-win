// Original C++ function from https://github.com/hms-dbmi/pagoda2
#include <RcppArmadillo.h>
#include <omp.h>
#include <progress.hpp>

using namespace std;
using namespace Rcpp;

// calculate column mean and variance, optionally taking a subset of rows to operate on
// [[Rcpp::export]]
Rcpp::DataFrame colMeanVarS(SEXP sY,  SEXP rowSel, int ncores=1) {
  // need to do this as SEXP, modify the slots on the fly
  S4 mat(sY);  
  const arma::uvec i(( unsigned int *)INTEGER(mat.slot("i")),LENGTH(mat.slot("i")),false,true); 
  const arma::ivec dims(INTEGER(mat.slot("Dim")),LENGTH(mat.slot("Dim")),false,true); 
  const arma::ivec p(INTEGER(mat.slot("p")),LENGTH(mat.slot("p")),false,true); 
  arma::vec Y(REAL(mat.slot("x")),LENGTH(mat.slot("x")),false,true); 
  
  bool rowSelSpecified=!Rf_isNull(rowSel);
  const arma::ivec rs=(rowSelSpecified) ? arma::ivec(INTEGER(rowSel),LENGTH(rowSel),false,true) : arma::ivec(); 
  
  int ncols=p.size()-1;
  int nrows=dims[0]; 
  if(rowSelSpecified) { 
    nrows=0; 
    for(int j=0;j<rs.size();j++) { if(rs[j]) { nrows++; } }
  }
  arma::vec meanV(ncols,arma::fill::zeros); arma::vec varV(ncols,arma::fill::zeros); arma::vec nobsV(ncols,arma::fill::zeros);
  // for each gene
#pragma omp parallel for num_threads(ncores) shared(meanV,varV,nobsV) 
  for(int g=0;g<ncols;g++) {
    int p0=p[g]; int p1=p[g+1]; 
    if(p1-p0 <1) { continue; }
    arma::colvec ly;
    if(rowSelSpecified) {
      // select valid rows
      int nvalid=0;
      ly=arma::vec(p1-p0);
      for(int j=p0;j<p1;j++) {
        if(rs[i[j]]) { 
          ly[nvalid]=Y[j]; nvalid++;
        }
      }
      nobsV[g]=nvalid;
      ly=ly.head(nvalid);
    } else {
      nobsV[g]=p1-p0;
      ly=Y.subvec(p0,p1-1);
    }
    
    double m=sum(ly)/nrows;
    meanV[g]=m;
    ly-=m; ly%=ly; 
    varV[g]=(sum(ly)+(m*m*(nrows-ly.size())))/nrows;
  }
  return Rcpp::DataFrame::create(Named("m")=meanV,Named("v")=varV,Named("nobs",nobsV));
}
