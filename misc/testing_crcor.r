library(mfutils)
library(data.table)
library(profvis)
library(here)

source(here("src", "stroop-rsa-pc.R"))

n_trial <- 100
n_vertex <- 100
n_resamps <- 10^(0:4)
expected_min <- 10
B1 <- matrix(rnorm(n_trial*n_vertex), nrow = n_trial, dimnames = list(NULL, rep(letters[1:10], 10)))
B2 <- matrix(rnorm(n_trial*n_vertex), nrow = n_trial, dimnames = list(NULL, rep(letters[1:10], 10)))


averaging_matrix <- function(x) {
    A <- indicator(x)
    As <- tcrossprod(A, diag(1/colSums(A)))
    colnames(As) <- unique(x)
    As
}


get_resampled_idx_new <- function(conditions, n_resamples, expected_min, seed = 0) {
    stopifnot(is.character(conditions) || is.numeric(n_resamples) || is.numeric(expected_min))
    set.seed(seed)

    n_conditions <- length(conditions)
    groups_list <- split(seq_along(conditions), conditions)
    resample_to <- Reduce(min, lapply(groups_list, length))
    if (resample_to != expected_min) stop("unexpected minimum n trials")
    
    idx <- t(replicate(n_resamples, unlist(lapply(groups_list, resample, size = resample_to))))
    colnames(idx) <- rep(names(groups_list), each = expected_min)
    idx
}

crcor1 <- function(x1, x2, n_resamples, expected_min){
    stopifnot(length(dim(x1)) == 2 || length(dim(x2)) == 2)
    
    resample_idx1 <- 
        get_resampled_idx_new(conditions = colnames(x1), n_resamples = n_resamples, expected_min = expected_min)
    resample_idx2 <- 
        get_resampled_idx_new(conditions = colnames(x2), n_resamples = n_resamples, expected_min = expected_min)

    A1 <- averaging_matrix(colnames(resample_idx1))
    A2 <- averaging_matrix(colnames(resample_idx2))

    res <- vector("list", n_resamples)
    for (ii in seq_len(n_resamples)) {
        idx1 <- resample_idx1[ii, ]
        idx2 <- resample_idx2[ii, ]
        x1i <- x1[, idx1, drop = FALSE]
        x2i <- x2[, idx2, drop = FALSE]
        res[[ii]] <- atanh(cor(x1i %*% A1, x2i %*% A2))
    }
    tanh(Reduce("+", res) / length(res))  ## take mean over resamples and invert atanh

}



get_resampled_idx_new2 <- function(conditions, n_resamples, expected_min, seed = 0) {
    stopifnot(is.character(conditions) || is.numeric(n_resamples) || is.numeric(expected_min))
    set.seed(seed)

    n_conditions <- length(conditions)
    groups_list <- split(seq_along(conditions), conditions)
    resample_to <- Reduce(min, lapply(groups_list, length))
    if (resample_to != expected_min) stop("unexpected minimum n trials")
    
    replicate(n_resamples, unlist(lapply(groups_list, resample, size = resample_to)), simplify = FALSE)

}


crcor2 <- function(x1, x2, n_resamples, expected_min){
    stopifnot(length(dim(x1)) == 2 || length(dim(x2)) == 2)
    
    nms1 <- colnames(x1)
    nms2 <- colnames(x2)
    g1 <- unique(nms1)
    g2 <- unique(nms2)

    resample_idx1 <- 
        get_resampled_idx_new2(conditions = nms1, n_resamples = n_resamples, expected_min = expected_min)
    resample_idx2 <- 
        get_resampled_idx_new2(conditions = nms2, n_resamples = n_resamples, expected_min = expected_min)

    A1 <- averaging_matrix(rep(g1, each = expected_min))
    A2 <- averaging_matrix(rep(g2, each = expected_min))
    
    res <- vector("list", n_resamples)
    for (ii in seq_len(n_resamples)) {
        idx1 <- resample_idx1[[ii]]
        idx2 <- resample_idx2[[ii]]
        x1i <- x1[, idx1, drop = FALSE]
        x2i <- x2[, idx2, drop = FALSE]
        res[[ii]] <- atanh(cor(x1i %*% A1, x2i %*% A2))
    }
    tanh(Reduce("+", res) / length(res))  ## take mean over resamples and invert atanh

}


crcor3 <- function(x1, x2, n_resamples, expected_min){
    stopifnot(length(dim(x1)) == 2 || length(dim(x2)) == 2)
    
    nms1 <- colnames(x1)
    nms2 <- colnames(x2)
    g1 <- unique(nms1)
    g2 <- unique(nms2)

    resample_idx1 <- 
        get_resampled_idx_new2(conditions = nms1, n_resamples = n_resamples, expected_min = expected_min)
    resample_idx2 <- 
        get_resampled_idx_new2(conditions = nms2, n_resamples = n_resamples, expected_min = expected_min)

    A1 <- averaging_matrix(rep(g1, each = expected_min))
    A2 <- averaging_matrix(rep(g2, each = expected_min))
    
    res <- matrix(0, ncol = ncol(A2), nrow = ncol(A1), dimnames = list(colnames(A1), colnames(A2)))
    for (ii in seq_len(n_resamples)) {
        idx1 <- resample_idx1[[ii]]
        idx2 <- resample_idx2[[ii]]
        x1i <- x1[, idx1, drop = FALSE]
        x2i <- x2[, idx2, drop = FALSE]
        res <- res + atanh(cor(x1i %*% A1, x2i %*% A2))
    }
    tanh(res / n_resamples)

}




lapply(n_resamps, function(x) system.time(crcor(B1, B2, x, expected_min = 10)))
lapply(n_resamps, function(x) system.time(crcor1(B1, B2, x, expected_min = 10)))
lapply(n_resamps, function(x) system.time(crcor2(B1, B2, x, expected_min = 10)))
lapply(n_resamps, function(x) system.time(crcor3(B1, B2, x, expected_min = 10)))

res <- crcor(B1, B2, 100, expected_min = 10)
res1 <- crcor1(B1, B2, 100, expected_min = 10)
res2 <- crcor2(B1, B2, 100, expected_min = 10)
res3 <- crcor3(B1, B2, 100, expected_min = 10)

identical(res, res1)
identical(res, res2)
identical(res, res3)

library(microbenchmark)
mbm = microbenchmark(
    res1 = crcor1(B1, B2, 10000, expected_min = 10),
    res2 = crcor2(B1, B2, 10000, expected_min = 10),
    res3 = crcor3(B1, B2, 10000, expected_min = 10),
    times=10
)
mbm
ggplot2::autoplot(mbm)
