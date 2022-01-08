# Function declaration
bootstrap_r <- function(ds, B = 1000) {
  # Preallocate storage for statistics
  boot_stat <- matrix(NA, nrow = B, ncol = 2)
  # Number of observations
  n <- length(ds)
  # Perform bootstrap
  for(i in seq_len(B)) {
    # Sample initial data
    gen_data <- ds[ sample(n, n, replace=TRUE) ]
    # Calculate sample data mean and SD
    boot_stat[i,] <- c(mean(gen_data),
                       sd(gen_data))
  }
  # Return bootstrap result
  return(boot_stat)
}

# Set seed to generate data
set.seed(512)
# Generate data
initdata <- rnorm(1000, mean = 21, sd = 10)
# Set a new _different_ seed for bootstrapping
set.seed(883)
# Perform bootstrap
result_r <- bootstrap_r(initdata)

Rcpp::sourceCpp("misc/bootstrap_cpp.cpp")

# Use the same seed use in R and C++
set.seed(883)
# Perform bootstrap with C++ function
result_cpp <- bootstrap_cpp(initdata)


all.equal(result_cpp, result_r)
norm(result_cpp - result_r)