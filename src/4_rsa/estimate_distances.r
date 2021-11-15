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
source(here("src", "stroop-rsa-pc.R"))

## set variables

task <- "Stroop"
n_resamples <- 1E3

if (interactive()) {  ## add variables (potentially unique to this script) useful for dev
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018Dev"
    prewh <- "none"
    measure <- "crcor"  ## crcor
    subjects <- fread(here("out/subjlist_ispc_retest.txt"))[[1]][1:5]
    waves <- c("wave1", "wave2")
    sessions <- "reactive"
    ttype_subset <- "bias"
    ii <- 1
    n_cores <- 10
    run_i <- 1
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

stopifnot(sessions %in% c("baseline", "proactive", "reactive"))
stopifnot(measure %in% c("crcor", "cveuc"))

if (sessions == "reactive" && ttype_subset == "bias") {
    expected_min <- 6
} else if (sessions == "proactive") {
    expected_min <- 3
}

atlas <- read_atlas(roiset)
rois <- names(atlas$roi)



## execute ----


input <- construct_filenames_h5(
    prefix = "coefs", subjects = subjects, waves = waves, sessions = sessions, rois = rois, runs = runs, 
    glmname = glmname, prewh = prewh
)
input[, g := paste0(subj, "__", wave, "__", session, "__", roi)]
l <- split(input, by = "g")

cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
res <- foreach(ii = seq_along(l), .final = function(x) setNames(x, names(l))) %dopar% {
    
    input_val <- l[[ii]]
    ses <- unique(input_val$session)

    B <- enlist(runs)
    for (run_i in seq_along(runs)) {
        Bi <- read_dset(input_val[run == runs[run_i], file_name], input_val[run == runs[run_i], dset_name])
        ttypes_to_get <- intersect(ttypes_by_run[[ses]][[run_i]], ttypes[[ttype_subset]])
        B[[run_i]] <- Bi[, colnames(Bi) %in% ttypes_to_get]  ## discard low-trialcount trialtypes
    }

    if (measure == "cveuc") {  ## cross-validated euclidean
        cvdist(
            average(B$run1, g = colnames(B$run1)),#[, ttypes_by_run[[ses]]$run1], ## paranoia, just to make sure all 
            average(B$run2, g = colnames(B$run2)),#[, ttypes_by_run[[ses]]$run2], ## conditions orders are the same
            scale = TRUE
        )
    } else if (measure == "crcor") {    ## cross-run correlation (with downsampling)
        crcor(B$run1, B$run2, n_resamples = n_resamples, expected_min = expected_min)
    }

}
stopImplicitCluster()


## write arrays within hdf5 file

fname <- construct_filename_rdm(measure = measure, glmname = glmname, roiset = roiset, prewh = prewh)

if (!file.exists(fname)) h5createFile(fname)
for (sub in subjects) {
    
    for (wav in waves) {

        for (ses in sessions) {
            
            nm <- paste0(sub, "__", wav, "__", ses)
            #group_name <- paste0("/", nm)            
            #if (!group_name %in% h5ls(fname)$group) h5createGroup(fname, group_name)

            ## extract data and concatenate into 3D array (RDM, roi)
            dat <- res[grep(nm, names(res))]
            dat <- abind(dat, rev.along = 0)
            dimnames(dat)[[3]] <- gsub(paste0(nm, "__"), "", dimnames(dat)[[3]])
            
            h5write(dat, fname, nm)

        }
    }
}

