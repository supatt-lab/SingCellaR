#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <vector>
#include <Rcpp.h>
#include "RcppPerpendicular.h"
#include "nn_parallel.h"
#include <kissrandom.h>
using namespace Rcpp;

#if defined(__MINGW32__)
#undef Realloc
#undef Free
#endif

#define __ERROR_PRINTER_OVERRIDE__ REprintf

//the compute_weights function was modified from the cytofkit package (Chen Hao, Date: 25/09/2015); 

// [[Rcpp::export]]
NumericVector compute_weights(NumericMatrix idx) {
  int nrow = idx.nrow(), ncol = idx.ncol();
  NumericMatrix weights(nrow*ncol, 3);
  int r = 0;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int k = idx(i,j)-1;
      weights(r, 0) = i+1;
      weights(r, 1) = k+1;
      NumericVector nodei = idx(i,_);
      NumericVector nodej = idx(k,_);
      double xi = intersect(nodei, nodej).size();  // count intersection number
      double yi= (2*ncol)-xi;  // count union number
      if(xi>0){
        double z = xi/yi;
        if(z > 0.05){
          weights(r, 2) = z;
          r++;
        }
      }
    }
  }
  return weights;
}

// [[Rcpp::export]]
NumericVector getDetectedGenesPerCell(arma::sp_mat X, int k) {
  int nc=X.n_cols;
  int nr=X.n_rows;
  NumericVector y(nc);
  for (int j = 0; j < nc;j++) {
    double sum = 0.0;
    for (int i = 0; i < nr;i++) {
      if(X(i,j) >= k){
        sum++;
      }
    }
    y[j] = sum;
  }
  return y;
}

// [[Rcpp::export]]
NumericVector getExpressingCellsPerGene(arma::sp_mat X, int k) {
  int nc=X.n_cols;
  int nr=X.n_rows;
  NumericVector y(nr);
  for (int i = 0; i < nr;i++) {
    double sum = 0.0;
    for (int j = 0; j < nc;j++) {
      if(X(i,j) >= k){
        sum++;
      }
    }
    y[i] = sum;
  }
  return y;
}

// [[Rcpp::export]]
arma::sp_mat makeBinaryMatrix(arma::sp_mat X,int k) {
  
  for(arma::sp_mat::const_iterator i = X.begin(); i != X.end(); ++i){
    if(X(i.row(),i.col()) >= k){
      X(i.row(),i.col()) = 1;
    }
  }
  return X;
}

// [[Rcpp::export]]
arma::sp_mat makeBinaryMatrix_for_float(arma::sp_mat X) {
  
  for(arma::sp_mat::const_iterator i = X.begin(); i != X.end(); ++i){
    if(X(i.row(),i.col()) > 0){
      X(i.row(),i.col()) = 1;
    }
  }
  return X;
}

// [[Rcpp::export]]
NumericVector count_expressing_cells(arma::sp_mat X,arma::sp_mat Y) {
  
  int Xcol=X.n_cols;
  int Xrow=X.n_rows;
  int Ycol=Y.n_cols;
  
  NumericMatrix t(Xrow, 2);
  for (int i = 0; i < Xrow;i++) {
    int sumX = 0;
    for (int j = 0; j < Xcol;j++) {
      if(X(i,j) > 0){
        sumX++;
      }
    }
    int sumY = 0;
    for (int k = 0; k < Ycol;k++) {
      if(Y(i,k) > 0){
        sumY++;
      }
    }
    t(i, 0) = sumX;
    t(i, 1) = sumY;
  }
  return t;
}

// [[Rcpp::export]]
arma::sp_mat count_divided_by_libsize(arma::sp_mat X) {
  
  int ncol = X.n_cols;
  int nrow = X.n_rows;
  arma::sp_mat Y(nrow,ncol);
  NumericVector z(ncol);
  for (int i = 0; i < ncol; ++i){
    z(i) = 0;
  }
  for(arma::sp_mat::const_iterator i = X.begin(); i != X.end(); ++i){
    z(i.col()) += *i;
  }
  for(arma::sp_mat::const_iterator i = X.begin(); i != X.end(); ++i){
    Y(i.row(),i.col()) = *i/z(i.col());
  }
  return Y;
}

//Below functions are copied from the uwot package.
//https://github.com/jlmelville/uwot/blob/master/src/nn_parallel.cpp

template <typename UwotAnnoyDistance>
auto annoy_nns_impl(const std::string &index_name, NumericMatrix mat,
                    std::size_t n_neighbors, std::size_t search_k,
                    std::size_t n_threads = 0, std::size_t grain_size = 1)
  -> List {
    
    std::size_t nrow = mat.rows();
    std::size_t ncol = mat.cols();
    
    std::vector<double> vmat = as<std::vector<double>>(mat);
    
    NNWorker<UwotAnnoyDistance> worker(index_name, vmat, ncol, n_neighbors,
                                       search_k);
    RcppPerpendicular::parallel_for(0, nrow, worker, n_threads, grain_size);
    
    return List::create(
      _("item") = IntegerMatrix(nrow, n_neighbors, worker.idx.begin()),
      _("distance") = NumericMatrix(nrow, n_neighbors, worker.dists.begin()));
  }

// [[Rcpp::export]]
List annoy_search_parallel_cpp(const std::string &index_name, NumericMatrix mat,
                               std::size_t n_neighbors, std::size_t search_k,
                               const std::string &metric,
                               std::size_t n_threads = 0,
                               std::size_t grain_size = 1) {
  if (metric == "euclidean") {
    return annoy_nns_impl<UwotAnnoyEuclidean>(index_name, mat, n_neighbors,
                                              search_k, n_threads, grain_size);
  } else if (metric == "cosine") {
    return annoy_nns_impl<UwotAnnoyCosine>(index_name, mat, n_neighbors,
                                           search_k, n_threads, grain_size);
  } else if (metric == "manhattan") {
    return annoy_nns_impl<UwotAnnoyManhattan>(index_name, mat, n_neighbors,
                                              search_k, n_threads, grain_size);
  } else if (metric == "hamming") {
    return annoy_nns_impl<UwotAnnoyHamming>(index_name, mat, n_neighbors,
                                            search_k, n_threads, grain_size);
  } else {
    stop("Unknown metric '", metric, "'");
  }
}