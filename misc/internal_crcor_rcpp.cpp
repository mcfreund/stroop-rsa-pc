#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat internal_crcor_rcpp(arma::mat x1,
               arma::mat x2,
               arma::umat idx1,
               arma::umat idx2,
               arma::mat A1,
               arma::mat A2,
               int n_resamples
               ) {
  // Preallocate storage for statistics
  int ncond1 = A1.n_cols;
  int ncond2 = A2.n_cols;
  arma::mat R(ncond1, ncond2);
  for(int i = 0; i < n_resamples; i++){
    // Get indices
    arma::uvec i1 = idx1.row(i).t();
    arma::uvec i2 = idx2.row(i).t();
    // Subset data
    arma::mat x1i = x1.cols(i1);
    arma::mat x2i = x2.cols(i2);
    // Average
    arma::mat x1i_bar = x1i * A1;
    arma::mat x2i_bar = x2i * A2;
    // Correlate
    R += arma::atanh(arma::cor(x1i_bar, x2i_bar));
  }
  return arma::tanh(R/n_resamples);
  }

