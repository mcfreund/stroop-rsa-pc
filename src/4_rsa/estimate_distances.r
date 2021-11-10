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

if (interactive()) { 
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018Dev"
    prewh <- "none"
    measure <- "cveuc"  ## crcor
    subjects <- c("115825", "130518", "155938", "178647")
    #subjects <- fread(here("out/subjlist_ispc_retest.txt"))[[1]]
    waves <- c("wave1", "wave2")
    sessions <- "proactive"
    ii <- 1
} else {
    args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
    glmname <- args$glmname
    roiset <- args$roiset
    prewh <- args$prewh
    measure <- args$measure
    outfile_prefix <- args$outfile_prefix
    subjects <- fread(here(args$subjlist))[[1]]
    waves <- strsplit(args$waves, " ")[[1]]
    sessions <- strsplit(args$sessions, " ")[[1]]
}

atlas <- read_atlas(roiset)
rois <- names(atlas$roi)
filename_master <- here("out", "parcellated", "master.h5")
info_master <- read_master(filename_master)[
    prefix == "coefs" & glm %in% paste0("glm-", glmname) & roiset == roiset & prewh %in% prewh &
    subj %in% subjects & wave %in% waves & session %in% sessions
    ]



## execute ----

info_master[, g := paste0(subj, "__", wave, "__", session, "__", roi)]
l <- split(info_master, by = "g")  ## NB: this calls split.data.table NOT split.data.frame! must use "by" arg.

cl <- makeCluster(n_core - 2, type = "FORK")
registerDoParallel(cl)
res <- foreach(ii = seq_along(l), .final = function(x) setNames(x, names(l))) %dopar% {
    
    input_val <- l[[ii]]
    sub <- unique(input_val$subj)
    wav <- unique(input_val$wave)
    ses <- unique(input_val$session)
    roi <- unique(input_val$roi)

    stopifnot(nrow(input_val) == 2)  ## one for each run
    
    B <- enlist(runs)
    for (run_i in seq_along(runs)) {
        Bi <- read_dset(filename_master, input_val[run == runs[run_i], name])
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


