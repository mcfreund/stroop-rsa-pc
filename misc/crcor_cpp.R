library(mfutils)
library(data.table)
library(here)
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(profvis)
source(here("src", "stroop-rsa-pc.R"))

n_trial <- 100
n_vertex <- 1000
n_resamples <- 10000
expected_min <- n_trial/10
set.seed(0)
x1 <- matrix(rnorm(n_trial*n_vertex), nrow = n_vertex, dimnames = list(NULL, rep(letters[1:10], n_trial/10)))
x2 <- matrix(rnorm(n_trial*n_vertex), nrow = n_vertex, dimnames = list(NULL, rep(letters[1:10], n_trial/10)))

get_resampled_idx <- function(conditions, n_resamples, expected_min, seed = 0) {
  stopifnot(is.character(conditions) || is.numeric(n_resamples) || is.numeric(expected_min))
  set.seed(seed)
  
  groups_list <- split(seq_along(conditions), conditions)
  resample_to <- Reduce(min, lapply(groups_list, length))
  if (resample_to != expected_min) stop("unexpected minimum n trials")
  
  t(replicate(n_resamples, unlist(lapply(groups_list, resample, size = resample_to))))
  
}


## functions ----

internal_crcor_r <- function(x1, x2, idx1, idx2, A1, A2, n_resamples) {
  res <- matrix(0, ncol = ncol(A2), nrow = ncol(A1), dimnames = list(colnames(A1), colnames(A2)))
  for (ii in seq_len(n_resamples)) {
    idx1_i <- idx1[ii, ]
    idx2_i <- idx2[ii, ]
    x1i <- x1[, idx1_i, drop = FALSE]
    x2i <- x2[, idx2_i, drop = FALSE]
    res <- res + atanh(cor(x1i %*% A1, x2i %*% A2))
  }
  tanh(res / n_resamples)
}

sourceCpp("misc/internal_crcor_rcpp.cpp")
sourceCpp("misc/internal_crcor_rcpp2.cpp")

crcor <- function(x1, x2, n_resamples, expected_min, f) {
  stopifnot(length(dim(x1)) == 2 || length(dim(x2)) == 2)
  
  nms1 <- colnames(x1)
  nms2 <- colnames(x2)
  g1 <- unique(nms1)
  g2 <- unique(nms2)
  
  idx1 <-
    get_resampled_idx(conditions = nms1, n_resamples = n_resamples, expected_min = expected_min)
  idx2 <-
    get_resampled_idx(conditions = nms2, n_resamples = n_resamples, expected_min = expected_min)
  
  A1 <- averaging_matrix(rep(g1, each = expected_min))
  A2 <- averaging_matrix(rep(g2, each = expected_min))
  
  if (f == "r") {
    res <- internal_crcor_r(x1, x2, idx1, idx2, A1, A2, n_resamples)
  } else if (f == "cpp") {
    res <- internal_crcor_rcpp(x1, x2, idx1 - 1L, idx2 - 1L, A1, A2, n_resamples)
  } else if (f == "cpp2") {
    res <- internal_crcor_rcpp2(x1, x2, idx1 - 1L, idx2 - 1L, A1, A2, n_resamples)
  }
  
  dimnames(res) <- list(g1, g2)
  
  res
  
}


## test ----

res_r <- crcor(x1, x2, n_resamples, expected_min, "r")
res_cpp <- crcor(x1, x2, n_resamples, expected_min, "cpp")
res_cpp2 <- crcor(x1, x2, n_resamples, expected_min, "cpp2")

all.equal(res_r, res_cpp)
all.equal(res_r, res_cpp2)


a <- microbenchmark(
  r     = crcor(x1, x2, n_resamples, expected_min, "r"),
  rcpp  = crcor(x1, x2, n_resamples, expected_min, "cpp"),
  rcpp2 = crcor(x1, x2, n_resamples, expected_min, "cpp2"),
  times = 20
)
ggplot2::autoplot(a)


p_r <- profvis(crcor(x1, x2, n_resamples, expected_min, "r"))
p_cpp <- profvis(crcor(x1, x2, n_resamples, expected_min, "cpp"))
p_cpp2 <- profvis(crcor(x1, x2, n_resamples, expected_min, "cpp2"))

htmlwidgets::saveWidget(p_r, "misc/profile_crcor_r.html")
htmlwidgets::saveWidget(p_cpp, "misc/profile_crcor_rcpp.html")
htmlwidgets::saveWidget(p_cpp2, "misc/profile_crcor_rcpp2.html")

