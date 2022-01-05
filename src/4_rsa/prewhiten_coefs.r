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
    atlas_name <- "schaefer2018_17_200"
    roi_col <- "parcel"
    space <- "fsaverage5"
    prewh <- "obsallave"
    subjlist <- "mi1"
    subjects <- fread(here(paste0("out/subjlist_", subjlist, ".txt")))[[1]]
    waves <- "wave1"
    sessions <- "proactive"
    n_cores <- 10
    ii <- 721
    overwrite <- FALSE
    n_resamples <- 1E2
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

stopifnot(prewh %in% expected$prewh)
if (prewh != "obsallave") stop("currently not configured for prewhitening other than 'obsallave'")

atlas <- load_atlas(atlas_name, space)
rois <- unique(atlas$key[[roi_col]])
## remove hippocampi from glasser atlas (not represented in fsaverages???)
if (atlas_name == "glasser2016" && grepl("fsaverage", space)) {
    rois <- rois[rois != "L_H" & rois != "R_H"]
}
roiset <- paste0(atlas_name, "_", roi_col)



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

input[, id := paste0(subj, wave, session, roi)]
id <- unique(input$id)

cl <- makeCluster(n_cores, type = "FORK", outfile = "")
registerDoParallel(cl)
res <- foreach(ii = seq_along(id), .inorder = FALSE) %dopar% {

    input_ii <- input[id == id[ii], ]  ## get all runs for single subj*wave*sess*roi
    
    subject <- input_ii[, unique(subj)]
    session <- input_ii[, unique(session)]
    wave <- input_ii[, unique(wave)]
    roi <- input_ii[, unique(roi)]

    ## read betas and extract good vertices (i.e., with signal variance), idx:
    B <- enlist(runs)
    idx <- enlist(runs)
    for (run in runs) {
        B[[run]] <- read_dset(input_ii$file_name[input_ii$run == run], input_ii$dset_name[input_ii$run == run])
        ## get inds for verts with BOLD:
        idx[[run]] <- which(!is_equal(Var(B[[run]], 1), 0))
    }
    idx <- Reduce(union, idx)
    B_good <- lapply(B, function(x) x[idx, ])

    ## remove trial-wise variance among conditions from each vertex:
    resids <- enlist(runs)
    for (run in runs) {
        X <- indicator(colnames(B_good[[run]]))  ## colnames of B gives the condition
        resids[[run]] <- resid(.lm.fit(X, t(B_good[[run]])))
    }
    
    ## estimate average covariance matrix:
    S <- lapply(resids, CovEst.2010OAS)
    S_bar <- Reduce("+", lapply(S, function(x) x$S)) / n_run  ## average over folds
    W <- sqrtm(S_bar)$Binv  ## invert and sqrt, yields sqrt of noise precision matrix; NB: sqrtm is slow in high-D

    ## apply to each run and save:
    out <- enlist(runs)
    for (run in runs) {    
        
        Bw <- crossprod(W, B_good[[run]])

        out[[run]] <- write_dset(
            mat = Bw,
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

    }
    
    for (run in runs) out[[run]] <- append(out[[run]], list(rho = S[[run]]$rho))
    do.call(rbind, out)  ## for each run, return file_name, dset_name, and shrinkage factor rho

}
stopCluster(cl)

data_info <- as.data.table(do.call(rbind, res))
fn <- construct_filename_datainfo(
        prefix = "coefs", subjlist = subjlist, glmname = glmname, roiset = roiset, prewh = prewh
        )
fwrite(data_info, fn)