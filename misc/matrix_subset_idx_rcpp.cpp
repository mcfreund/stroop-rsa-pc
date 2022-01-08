#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_subset_idx_rcpp(
    Rcpp::NumericMatrix x, Rcpp::IntegerVector y) { 
  
  // Determine the number of observations
  int n_cols_out = y.size();
  
  // Create an output matrix
  Rcpp::NumericMatrix out = Rcpp::no_init(x.nrow(), n_cols_out);
  
  // Loop through each column and copy the data. 
  for(unsigned int z = 0; z < n_cols_out; ++z) {
    out(Rcpp::_, z) = x(Rcpp::_, y[z]);
  }
  
  return out;
}