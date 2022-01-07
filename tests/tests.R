
# library(colorout)
# library(here)
# library(dplyr)
# library(tidyr)
# library(data.table)
# library(gifti)
# library(abind)



library(colorout)
library(here)
source(here("src", "stroop-rsa-pc.R"))

atlas <- read_atlas()
files <- c(
    "/data/nil-external/ccp/freund/stroop-rsa-pc/out/glms/115825/wave1/RESULTS/Stroop/reactive_lsall_1rpm/STATS_115825_run1_L_REML.func.gii",
    "/data/nil-external/ccp/freund/stroop-rsa-pc/out/glms/115825/wave1/RESULTS/Stroop/reactive_lsall_1rpm/STATS_115825_run1_R_REML.func.gii",
    "/data/nil-external/ccp/freund/stroop-rsa-pc/out/glms/115825/wave1/RESULTS/Stroop/reactive_lsall_1rpm/STATS_115825_run2_L_REML.func.gii",
    "/data/nil-external/ccp/freund/stroop-rsa-pc/out/glms/115825/wave1/RESULTS/Stroop/reactive_lsall_1rpm/STATS_115825_run2_R_REML.func.gii"
)
l <- lapply(files, gifti::read_gifti)
gii <- l[[1]]

new_parcellated_data()
str(extract_labels(gii))
str(extract_data(gii, pattern = "_Coef"))
str(concat_hemis(l[1:2], pattern = "_Coef"))
str(parcellate_data(concat_hemis(l[1:2], pattern = "_Coef"), atlas))
str(parcellate_data(concat_hemis(l[1:2], pattern = "_Coef"), atlas))
x <- giftis_to_parcellated_image(
    giftis = l, folds = c("run1", "run1", "run2", "run2"), hemis = c("L", "R", "L" , "R"), atlas = atlas,
    glm_name = "lsall_1rpm", roi_set = "Schaefer2018_control", prewhitened = "none", 
    shrinkage_var = numeric(),
    shrinkage_cov = numeric(),
    subject = "",
    wave = "",
    session = "",
    task = "",
    pattern = "_Coef"
)
x
lapply(x$data, class)
x$bad_vertices


file.exists(construct_filename_gifti("132017", "wave1", "baseline", "run1", "lsall_1rpm", "L"))
file.exists(construct_filename_gifti("132017", "wave1", c("baseline", "baseline"), "run1", "lsall_1rpm", "L"))

file.exists(construct_filename_h5("pimage-B", "lsall_1rpm", "Schaefer2018_control", "none"))
file.exists(construct_filename_h5("pimage-B", "lsall_1rpm", "Schaefer2018_control", c("none", "none")))

create_nested_group()


a <- matrix(rnorm(1E6), ncol = 200)
nms <- write_dset(
    a, dset_prefix = "coefs",
    subject = "132017", wave = "wave1", session = "baseline", run = "run1", 
    roiset = "Schaefer2018_control", roi = "Vis",
    glmname = "lsall-1rpm", prewh = "none"
    )
write_links(list(nms))
h5ls(here::here("out", "parcellated", "master.h5"))


## -----
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

## original:

get_resampled_idx <- function(conditions, n_resamples, expected_min, seed = 0) {
    stopifnot(is.character(conditions) || is.numeric(n_resamples) || is.numeric(expected_min))
    set.seed(seed)

    n_conditions <- length(conditions)
    groups_list <- split(seq_along(conditions), conditions)
    resample_to <- Reduce(min, lapply(groups_list, length))
    if (resample_to != expected_min) stop("unexpected minimum n trials")
    
    t(replicate(n_resamples, unlist(lapply(groups_list, resample, size = resample_to))))

}
crcor <- function(x1, x2, n_resamples, expected_min){
    stopifnot(length(dim(x1)) == 2 || length(dim(x2)) == 2)
    
    resample_idx1 <- 
        get_resampled_idx(conditions = colnames(x1), n_resamples = n_resamples, expected_min = expected_min)
    resample_idx2 <- 
        get_resampled_idx(conditions = colnames(x2), n_resamples = n_resamples, expected_min = expected_min)

    n_resamples <- nrow(resample_idx1)
    res <- vector("list", n_resamples)
    for (ii in seq_len(n_resamples)) {
        idx1 <- resample_idx1[ii, ]
        idx2 <- resample_idx2[ii, ]
        x1i <- x1[, idx1, drop = FALSE]
        x2i <- x2[, idx2, drop = FALSE]
        res[[ii]] <- atanh(cor(average(x1i, colnames(x1i)), average(x2i, colnames(x2i))))
    }
    
    tanh(Reduce("+", res) / length(res))  ## take mean over resamples and invert atanh

}