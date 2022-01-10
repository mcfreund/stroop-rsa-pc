library(here)
library(data.table)
library(mfutils)
source(here("src", "stroop-rsa-pc.R"))

n_resamples <- 100

get_resampled_idx_orig <- function(conditions, n_resamples, expected_min, seed = 0) {
    stopifnot(is.character(conditions) || is.numeric(n_resamples) || is.numeric(expected_min))
    set.seed(seed)

    n_conditions <- length(conditions)
    groups_list <- split(seq_along(conditions), conditions)
    resample_to <- Reduce(min, lapply(groups_list, length))
    if (resample_to != expected_min) stop("unexpected minimum n trials")
    
    t(replicate(n_resamples, unlist(lapply(groups_list, resample, size = resample_to))))

}
crcor_orig <- function(x1, x2, n_resamples, expected_min){
    stopifnot(length(dim(x1)) == 2 || length(dim(x2)) == 2)
    
    resample_idx1 <- 
        get_resampled_idx_orig(conditions = colnames(x1), n_resamples = n_resamples, expected_min = expected_min)
    resample_idx2 <- 
        get_resampled_idx_orig(conditions = colnames(x2), n_resamples = n_resamples, expected_min = expected_min)

    n_resamples <- nrow(resample_idx1)
    res <- vector("list", n_resamples)
    for (ii in seq_len(n_resamples)) {
        idx1 <- resample_idx1[ii, ]
        idx2 <- resample_idx2[ii, ]
        x1i <- x1[, idx1, drop = FALSE]
        x2i <- x2[, idx2, drop = FALSE]
        res[[ii]] <- atanh(cor(mfutils::average(x1i, colnames(x1i)), mfutils::average(x2i, colnames(x2i))))
    }
    
    tanh(Reduce("+", res) / length(res))  ## take mean over resamples and invert atanh

}


fname1 <- here("tests", "DMCC2834766_wave1_proactive_run1_glasser2016_parcel_L_V8.h5")
fname2 <- here("tests", "DMCC2834766_wave1_proactive_run2_glasser2016_parcel_L_V8.h5")
dname1 <- "coefs__DMCC2834766__wave1__proactive__run1__glasser2016_parcel__L_V8__glm-lsall_1rpm__prewh-none"
dname2 <- "coefs__DMCC2834766__wave1__proactive__run2__glasser2016_parcel__L_V8__glm-lsall_1rpm__prewh-none"

B1 <- read_dset(fname1, dname1)
B2 <- read_dset(fname2, dname2)

res_orig <- crcor_orig(B1, B2, n_resamples, 3)
idx1 <- get_resampled_idx(colnames(B1), n_resamples)
idx2 <- get_resampled_idx(colnames(B2), n_resamples)
res <- crcor(x1 = B1, x2 = B2, idx1, idx2)
all.equal(res_orig, res)
#res_sort <- res[rownames(res_orig), colnames(res_orig)]
# colnames(res_orig) == colnames(res)
# rownames(res_orig) == rownames(res)
#norm(res_orig - res)
