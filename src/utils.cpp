#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppParallel)]]
//[[Rcpp::depends(RcppAnnoy)]]
#include <RcppParallel.h>
#include <annoylib.h>
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

template<typename S, typename T, typename Distance, typename Random>
struct NNWorker : public RcppParallel::Worker {
  std::string index_name;
  RcppParallel::RMatrix<double> mat;
  RcppParallel::RMatrix<double> dists;
  RcppParallel::RMatrix<int> idx;
  std::size_t ncol;
  std::size_t n_neighbors;
  std::size_t search_k;
  
  NNWorker(
    const std::string& index_name,
    const Rcpp::NumericMatrix& mat,
    Rcpp::NumericMatrix& dists,
    Rcpp::IntegerMatrix& idx,
    std::size_t ncol,
    std::size_t n_neighbors,
    std::size_t search_k
  ) :
    index_name(index_name), mat(mat), dists(dists), idx(idx), ncol(ncol), 
    n_neighbors(n_neighbors), search_k(search_k)
  {}
  
  void operator()(std::size_t begin, std::size_t end) {
    AnnoyIndex<S, T, Distance, Random> index(ncol);
    index.load(index_name.c_str());
    
    for (std::size_t i = begin; i < end; i++) {
      RcppParallel::RMatrix<double>::Row row = mat.row(i);
      std::vector<T> fv(row.length());
      std::copy(row.begin(), row.end(), fv.begin());
      std::vector<S> result;
      std::vector<T> distances;
      
      index.get_nns_by_vector(fv.data(), n_neighbors, search_k, &result, &distances);
      if (result.size() != n_neighbors || distances.size() != n_neighbors) { 
        break;
      }
      
      for (std::size_t j = 0; j < n_neighbors; j++) {
        dists(i, j) = distances[j];
        idx(i, j) = result[j];
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::List annoy_euclidean_nns(const std::string& index_name,
                               const Rcpp::NumericMatrix& mat,
                               std::size_t n_neighbors, std::size_t search_k,
                               std::size_t grain_size = 1,
                               bool verbose = false) {
  std::size_t nrow = mat.rows();
  std::size_t ncol = mat.cols();
  Rcpp::NumericMatrix dist(nrow, n_neighbors);
  Rcpp::IntegerMatrix idx(nrow, n_neighbors);
  idx.fill(-1);
  
  NNWorker<int32_t, float, Euclidean, Kiss64Random>
    worker(index_name, mat, dist, idx, ncol, n_neighbors, search_k);
  RcppParallel::parallelFor(0, nrow, worker, grain_size);
  
  return Rcpp::List::create(Rcpp::Named("item") = idx,
                            Rcpp::Named("distance") = dist);
}

// [[Rcpp::export]]
Rcpp::List annoy_cosine_nns(const std::string& index_name,
                            const Rcpp::NumericMatrix& mat,
                            std::size_t n_neighbors, std::size_t search_k,
                            std::size_t grain_size = 1,
                            bool verbose = false) {
  std::size_t nrow = mat.rows();
  std::size_t ncol = mat.cols();
  Rcpp::NumericMatrix dist(nrow, n_neighbors);
  Rcpp::IntegerMatrix idx(nrow, n_neighbors);
  idx.fill(-1);
  
  NNWorker<int32_t, float, Angular, Kiss64Random>
    worker(index_name, mat, dist, idx, ncol, n_neighbors, search_k);
  RcppParallel::parallelFor(0, nrow, worker, grain_size);
  
  return Rcpp::List::create(Rcpp::Named("item") = idx,
                            Rcpp::Named("distance") = dist);
}