#!/usr/bin/env bash Rscript --vanilla

## --------------------------------------------------------------------------------------------------------------------
## Script: 
## Author: Mike Freund (mfreundc@gmail.com)
## Input: 
##  - ...
## Output:
##  - ...
## Notes:
##  - ...
## --------------------------------------------------------------------------------------------------------------------


## setup ----


## packages and sourced variables

library(colorout)
library(here)
library(dplyr)
library(tidyr)
library(data.table)
library(gifti)
library(abind)
library(rhdf5)
library(foreach)
library(doParallel)
library(CovTools)
library(pracma)
source(here("src", "stroop-rsa-pc.R"))

## set variables

task <- "Stroop"
n_resample <- 1E4

if (interactive()) { 
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018Dev"
    prewh <- "obsall"  ## obsresamp, obsall
    subjects <- c("115825", "130518", "155938", "178647")
    subjlist <- "ispc_retest"
    subjects <- fread(here("out/subjlist_ispc_retest.txt"))[[1]]
    waves <- c("wave1", "wave2")
    sessions <- "reactive"
    ii <- 1
} else {
    args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
    glmname <- args$glmname
    roiset <- args$roiset
    prewh <- args$prewh
    subjects <- fread(here("out", paste0("subjlist_", args$subjlist, ".txt")))[[1]]
    waves <- strsplit(args$waves, " ")[[1]]
    sessions <- strsplit(args$sessions, " ")[[1]]
    n_cores <- args$n_cores
}

atlas <- read_atlas(roiset)
rois <- names(atlas$roi)



## execute ----

input <- construct_filenames_h5(
    prefix = "coefs", subjects = subjects, waves = waves, sessions = sessions, rois = rois, runs = runs, 
    glmname = glmname, prewh = prewh
)

cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
res <- foreach(ii = seq_along(input$file_name), .inorder = FALSE) %dopar% {
    
    input_val <- input[ii, ]
    B <- read_dset(input_val$file_name, input_val$dset_name)
    
    ## remove condition variance from trial-level obs (subtract condition mean pattern from each trial's pattern)
    X <- model.matrix(~ colnames(B) + 0)  ## colnames of B gives the condition
    resids <- resid(.lm.fit(X, B))
    ## estimate covariance matrix of residuals and invert
    if (prewh == "obsall") {
        S <- CovEst.2010OAS(t(resids))$S
    } else if (prewh == "obsresamp") {
        S <- resample_apply_combine(
            x = resids, 
            resample_idx = get_resampled_idx(conditions = colnames(resids), n_resample, expected_min = 1),
            apply_fun = function(.x) CovEst.2010OAS(t(.x))$S
            )
    }
    W <- sqrtm(S)$Binv  ## sqrt of inverse
    B_white <- crossprod(W, B)
    
    write_dset(
        mat = B_white, 
        dset_prefix = "coefs", 
        subject = input_val[, subj], 
        session = input_val[, session], 
        wave = input_val[, wave], 
        run = input_val[, run], 
        roiset = roiset, 
        roi = input_val[, roi],
        glmname = glmname,
        prewh = prewh, 
        write_colnames = TRUE
        )

}
stopImplicitCluster()

data_info <- as.data.table(do.call(rbind, res))
fn <- construct_filename_datainfo(
        prefix = "coefs", subjlist = subjlist, glmname = glmname, roiset = roiset, prewh = prewh
        )
fwrite(data_info, fn)