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
suppressMessages(library(here))
suppressMessages(library(dplyr))
library(tidyr)
suppressMessages(library(data.table))
suppressMessages(library(gifti))
library(abind)
library(rhdf5)
library(foreach)
suppressMessages(library(doParallel))
source(here("src", "stroop-rsa-pc.R"))

## set variables

task <- "Stroop"

if (interactive()) {  ## add variables (potentially unique to this script) useful for dev
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018Network"
    prewh <- "none"
    measure <- "crcor"  ## crcor
    subjects <- fread(here("out/subjlist_mimc.txt"))[[1]]
    waves <- "wave1"
    sessions <- "proactive"
    ttype_subset <- "pc50"
    ii <- 1
    n_cores <- 1
    run_i <- 1
    n_resamples <- 1
    overwrite <- FALSE
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

stopifnot(sessions %in% c("baseline", "proactive", "reactive"))
stopifnot(measure %in% c("crcor", "cveuc"))

atlas <- read_atlas(roiset)
rois <- names(atlas$roi)



## execute ----


input <- construct_filenames_h5(
    prefix = "coefs", subjects = subjects, waves = waves, sessions = sessions, rois = rois, runs = runs, 
    glmname = glmname, prewh = prewh
)
input[, g := paste0(subj, "__", wave, "__", session, "__", roi)]
l <- split(input, by = "g")

## skip computation of existing output files if overwrite == FALSE:

fname <- construct_filename_rdm(
    measure = measure, glmname = glmname, ttype_subset = ttype_subset, roiset = roiset, prewh = prewh
    )  ## output file

if (!overwrite) {
    if (file.exists(fname)) {
        group_name <- combo_paste(subjects, waves, sessions, sep = "__")
        group_exists <- group_name %in% h5ls(fname)$name
        if (any(group_exists)) {
            existing_group_idx <- grep(paste0(group_name[group_exists], collapse = "|"), names(l))
            l <- l[-existing_group_idx]
            if (length(l) == 0L) stop("All to-be-written files already exist. Exiting...")
        }
    }
}




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
            average(B$run1, g = colnames(B$run1)),
            average(B$run2, g = colnames(B$run2)),
            scale = TRUE
        )
    } else if (measure == "crcor") {    ## cross-run correlation (with downsampling)
        crcor(
            B$run1, B$run2, 
            n_resamples = n_resamples, 
            expected_min = expected_min[[paste0(ses, "_", ttype_subset)]]
            )
    }

}
stopCluster(cl)


## write arrays within hdf5 file

if (!file.exists(fname)) h5createFile(fname)
outnames <- unique(rbindlist(l)[, .(subj, wave, session)])
for (row_i in seq_len(nrow(outnames))) {
            
    nm <- paste0(outnames[row_i, ], collapse = "__")

            ## extract data and concatenate into 3D array (RDM, roi)
            dat <- res[grep(nm, names(res))]
            dat <- abind(dat, rev.along = 0)
            dimnames(dat)[[3]] <- gsub(paste0(nm, "__"), "", dimnames(dat)[[3]])
            
            h5write(dat, fname, nm)

        }