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
library(mfutils)
source(here("src", "stroop-rsa-pc.R"))
#library(profvis)

## set variables

task <- "Stroop"

if (interactive()) { 
    glmname <- "lsall_1rpm"
    atlas_name <- "glasser2016"
    roi_col <- "parcel"
    space <- "fsaverage5"
    prewh <- "obsall"  ## obsresamp, obsall, obsresampbias, obsresamppc50, obsbias, obspc50
    subjlist <- "mi1"
    subjects <- fread(here(paste0("out/subjlist_", subjlist, ".txt")))[[1]]
    waves <- "wave1"
    sessions <- "proactive"
    n_cores <- 10
    ii <- 2#321  ## Vis: 331, SomMot: 341
    overwrite <- FALSE
    n_resamples <- 1E2
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

stopifnot(prewh %in% expected$prewh)

atlas <- load_atlas(atlas_name, space)
rois <- unique(atlas$key[[roi_col]])
## remove hippocampi from glasser atlas (not represented in fsaverages???)
if (atlas_name == "glasser2016" && grepl("fsaverage", space)) {
    rois <- rois[rois != "L_H" & rois != "R_H"]
}
roiset <- paste0(atlas_name, "_", roi_col)


if (prewh %in% c("obsresamp", "obsall")) {
    ttype_subset <- "all"
} else if (prewh %in% c("obsresampbias", "obsbias")) {
    ttype_subset <- "bias"
} else if (prewh %in% c("obsresamppc50", "obspc50")) {
    ttype_subset <- "pc50"
} else {
    stop("unexpected input for prewh argument")
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
    subject <- input_val[, subj]
    session <- input_val[, session]
    wave <- input_val[, wave]
    run <- input_val[, run]
    roi <- input_val[, roi]

    B <- read_dset(input_val$file_name, input_val$dset_name)
    
    ## remove trial-wise variance due to conditions from each voxel
    X <- indicator(colnames(B))  ## colnames of B gives the condition
    resids <- resid(.lm.fit(X, t(B)))
    
    ## estimate covariance matrix of residuals and invert
    resids <- resids[rownames(resids) %in% ttypes[[ttype_subset]], ]  ## extract specified trialtypes
    if (grepl("resamp", prewh)) {
        S <- resample_apply_combine(
            x = t(resids), 
            resample_idx = get_resampled_idx(
                conditions = rownames(resids), 
                n_resamples, 
                expected_min = expected_min[[paste0(session, "_", ttype_subset)]]
                ),
            apply_fun = function(.x) CovEst.2010OAS(t(.x))$S,
            combine_fun = "iterative_add",
            outdim = c(ncol(resids), ncol(resids))
            )
    } else {
        S <- CovEst.2010OAS(resids)
    }
    
    W <- crossprod(sqrtm(S$S)$Binv, B)  ## apply sqrt of inverse

    out <- write_dset(
        mat = W,
        dset_prefix = "coefs", 
        subject = subject, 
        session = session, 
        wave = wave, 
        run = run, 
        roiset = roiset, 
        roi = roi,
        glmname = glmname,
        prewh = prewh, 
        write_colnames = TRUE
        )

    c(out, rho = S$rho)  ## returns metadata


}
stopCluster(cl)

data_info <- as.data.table(do.call(rbind, res))
fn <- construct_filename_datainfo(
        prefix = "coefs", subjlist = subjlist, glmname = glmname, roiset = roiset, prewh = prewh
        )
fwrite(data_info, fn)