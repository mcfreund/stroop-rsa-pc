#!/usr/bin/env bash Rscript --vanilla

## --------------------------------------------------------------------------------------------------------------------
## Script: 
## Author: Mike Freund (mfreundc@gmail.com)
## Input: 
##  - ...
## Output:
##  - ...
## Notes:
##  - TODO: should probably add handling of vertices with no bold
## --------------------------------------------------------------------------------------------------------------------


## setup ----


## packages and sourced variables

library(httpgd)
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
hgd()

## set variables

task <- "Stroop"
n_resample <- 1E4

if (interactive()) { 
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018Dev"
    prewh <- "obsall"  ## obsresamp, obsall
    subjlist <- "ispc_retest"
    subjects <- fread(here("out/subjlist_ispc_retest.txt"))[[1]][1:5]
    waves <- c("wave1", "wave2")
    sessions <- "reactive"
    n_cores <- 10
    ii <- 2#321  ## Vis: 331, SomMot: 341
    expected_min <- 1
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

atlas <- read_atlas(roiset)
rois <- names(atlas$roi)



## execute ----

input <- construct_filenames_h5(
    prefix = "coefs", subjects = subjects, waves = waves, sessions = sessions, rois = rois, runs = runs, 
    glmname = glmname, prewh = "none"
)

cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
res <- foreach(ii = seq_along(input$file_name), .inorder = FALSE) %dopar% {
    
    input_val <- input[ii, ]
    B <- read_dset(input_val$file_name, input_val$dset_name)
    
    ## remove trial-wise variance due to conditions from each voxel
    X <- indicator_matrix(colnames(B))  ## colnames of B gives the condition
    resids <- resid(.lm.fit(X, t(B)))
    
    ## estimate covariance matrix of residuals and invert
    if (prewh == "obsall") {
        S <- CovEst.2010OAS(resids)$S
    } else if (prewh == "obsresamp") {
        S <- resample_apply_combine(
            x = t(resids), 
            resample_idx = get_resampled_idx(conditions = rownames(resids), n_resample, expected_min = expected_min),
            apply_fun = function(.x) CovEst.2010OAS(t(.x))$rho
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