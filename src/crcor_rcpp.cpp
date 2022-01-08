#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat crcor_rcpp(arma::mat x1,
                              arma::mat x2,
                              arma::umat idx1,
                              arma::umat idx2,
                              arma::mat A1,
                              arma::mat A2
) {
  // Preallocate storage for statistics
  int ncond1 = A1.n_cols;
  int ncond2 = A2.n_cols;
  int n_resamples = idx1.n_rows;
  arma::mat R(ncond1, ncond2);
  for(int i = 0; i < n_resamples; i++){
    // Subset data, average, and correlate
    R += arma::atanh(arma::cor(x1.cols(idx1.row(i).t()) * A1, x2.cols(idx2.row(i).t()) * A2));
  }
  return arma::tanh(R/n_resamples);
}

