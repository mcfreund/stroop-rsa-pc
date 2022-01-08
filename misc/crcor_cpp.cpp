#include <Rcpp.h>

// https://stackoverflow.com/questions/62118084/rcpp-select-subset-numericmatrix-column-by-a-numericvector
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


// Function declaration with export tag
// [[Rcpp::export]]
Rcpp::NumericMatrix
  internal_crcor_cpp(Rcpp::NumericMatrix x1,
                Rcpp::NumericMatrix x2,
                Rcpp::NumericMatrix A1,
                Rcpp::NumericMatrix A2,
                Rcpp::IntegerMatrix idx1,
                Rcpp::IntegerMatrix idx2,
                int n_resamples = 10000) {
    int ncond1 = A1.ncol();
    int ncond2 = A2.ncol();
    // Preallocate storage for statistics
    Rcpp::NumericMatrix results(ncond1, ncond2);
    
    
    // Subset data matrices
    Rcpp::NumericMatrix x1_i = matrix_subset_idx_rcpp(x1, idx1( 0, Rcpp::_));
    Rcpp::NumericMatrix x2_i = matrix_subset_idx_rcpp(x2, idx2( 0, Rcpp::_));
    // Rcpp::NumericMatrix boot_stat(B, 2);
    // // Number of observations
    // int n = ds.size();
    // // Perform bootstrap
    // for(int i = 0; i < B; i++) {
    //   // Sample initial data
    //   Rcpp::NumericVector gen_data =
    //     ds[ floor(Rcpp::runif(n, 0, n)) ];
    //   // Calculate sample mean and std dev
    //   boot_stat(i, 0) = mean(gen_data);
    //   boot_stat(i, 1) = sd(gen_data);
    // }
    // Return bootstrap results
    return x1_i;
  }






// internal_crcor_r <- function(x1, x2, A1, A2, idx1, idx2, n_resamples) {
//   res <- matrix(0, ncol = ncol(A2), nrow = ncol(A1), dimnames = list(colnames(A1), colnames(A2)))
//   for (ii in seq_len(n_resamples)) {
//     idx1 <- resample_idx1[[ii]]
//     idx2 <- resample_idx2[[ii]]
//     x1i <- x1[, idx1, drop = FALSE]
//     x2i <- x2[, idx2, drop = FALSE]
//     res <- res + atanh(cor(x1i %*% A1, x2i %*% A2))
//   }
//   tanh(res / n_resamples)
// }
