#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <progress.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]


// [[Rcpp::export]]
arma::mat armaColSumFull(const arma::mat& m, int ncores=1, bool verbose=false)
{
  arma::vec s(m.n_cols, arma::fill::zeros);
  int tot = m.n_cols;
  Progress p(tot, verbose);
#pragma omp parallel for num_threads(ncores) shared(s)
  for(int i=0;i<(m.n_cols);i++) {
    s(i) = arma::sum(m.col(i));
  }
  return(s);
}

// [[Rcpp::export]]
arma::mat armaColSumSparse(const arma::sp_mat& m, int ncores=1, bool verbose=false)
{
  arma::vec s(m.n_cols, arma::fill::zeros);
  int tot = m.n_cols;
  Progress p(tot, verbose);
#pragma omp parallel for num_threads(ncores) shared(s)
  for(int i=0;i<(m.n_cols);i++) {
    s(i) = arma::sum(m.col(i));
  }
  return(s);
}

// calculate column mean and variance
// [[Rcpp::export]]
Rcpp::DataFrame colMeanVarS(const arma::sp_mat& m, int ncores=1) {
  
  arma::vec meanV(m.n_cols,arma::fill::zeros); arma::vec varV(m.n_cols,arma::fill::zeros); arma::vec nobsV(m.n_cols,arma::fill::zeros);
#pragma omp parallel for num_threads(ncores) shared(meanV,varV,nobsV)
  for(int i=0;i<(m.n_cols);i++) {
    meanV(i) = arma::mean(m.col(i));
    nobsV(i) = m.col(i).n_nonzero;
    varV(i) = arma::var(m.col(i));
  }
  return Rcpp::DataFrame::create(Named("m")=meanV,Named("v")=varV,Named("nobs",nobsV));
}
