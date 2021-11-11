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

if (interactive()) {  ## add variables (potentially unique to this script) useful for dev
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018Dev"
    prewh <- "none"
    measure <- "cveuc"  ## crcor
    subjects <- fread(here("out/subjlist_ispc_retest.txt"))[[1]][1:5]
    waves <- c("wave1", "wave2")
    sessions <- "reactive"
    ii <- 1
    n_cores <- 10
    run_i <- 1
} else {
    source(here("src", "parse_args.r"))
    print(args)
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
        Bi_bar <- average(Bi, g = colnames(Bi))
        colorder <- ttypes_by_run[[ses]][[run_i]]  ## paranoia, just to make sure all conditions orders are the same.
        B[[run_i]] <- Bi_bar[, colorder]
    }
    
    if (measure == "cveuc") {
        
        cvdist(B[[1]], B[[2]])  ## cross-validated euclidean distance

    } else if (measure == "crcor") {
        
        crcor(B[[1]], B[[2]])  ## downsampled cross-run correlation
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

