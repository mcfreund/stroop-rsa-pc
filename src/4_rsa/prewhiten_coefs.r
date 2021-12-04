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

suppressMessages(library(here))
suppressMessages(library(dplyr))
library(tidyr)
suppressMessages(library(data.table))
suppressMessages(library(gifti))
library(abind)
library(rhdf5)
library(foreach)
suppressMessages(library(doParallel))
suppressMessages(library(CovTools))
library(pracma)
source(here("src", "stroop-rsa-pc.R"))

## set variables

task <- "Stroop"

if (interactive()) { 
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018Dev"
    prewh <- "obsresampbias"  ## obsresamp, obsall, obsresampbias, obsresamppc50
    subjlist <- "ispc_retest"
    subjects <- fread(here("out/subjlist_ispc_retest.txt"))[[1]][1:5]
    waves <- c("wave1", "wave2")
    sessions <- "reactive"
    n_cores <- 10
    ii <- 2#321  ## Vis: 331, SomMot: 341
    overwrite <- FALSE
    n_resamples <- 1E2
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

stopifnot(prewh %in% expected$prewh)

atlas <- read_atlas(roiset)
rois <- names(atlas$roi)

if (prewh == "obsresamp") {
    ttype_subset <- "all"
} else if (prewh == "obsresampbias") {
    ttype_subset <- "bias"
} else if (prewh == "obsresamppc50") {
    ttype_subset <- "pc50"
}


## execute ----

input <- construct_filenames_h5(
    prefix = "coefs", subjects = subjects, waves = waves, sessions = sessions, rois = rois, runs = runs, 
    glmname = glmname, prewh = "none"
)

## option to skip if dset already exists (not overwrite):
if (!overwrite) {
    dset_name_new <- gsub("prewh-none", paste0("prewh-", prewh), input$dset_name)
    file_ls <- mclapply(input$file_name, h5ls, mc.cores = 20)
    dset_exists <- mapply(function(x, y) x %in% y$name, x = dset_name_new, y = file_ls)
    input <- input[!dset_exists, ]
}


cl <- makeCluster(n_cores, type = "FORK", outfile = "")
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

    } else if (grepl("resamp", prewh)) {

        resids <- resids[rownames(resids) %in% ttypes[[ttype_subset]], ]

        S <- resample_apply_combine(
            x = t(resids), 
            resample_idx = get_resampled_idx(
                conditions = rownames(resids), 
                n_resamples, 
                expected_min = expected_min[[paste0(ses, "_", ttype_subset)]])
                ),
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
stopCluster(cl)

data_info <- as.data.table(do.call(rbind, res))
fn <- construct_filename_datainfo(
        prefix = "coefs", subjlist = subjlist, glmname = glmname, roiset = roiset, prewh = prewh
        )
fwrite(data_info, fn)