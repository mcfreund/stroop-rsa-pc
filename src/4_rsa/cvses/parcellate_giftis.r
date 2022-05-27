#!/usr/bin/env bash Rscript --vanilla

## --------------------------------------------------------------------------------------------------------------------
## Script: 
## Author: Mike Freund (mfreundc@gmail.com)
## Input: 
##  - character vectors: glmname, roiset, subjlist, waves, sessions
##      - length one: glmname, roiset
## Output:
##  - vertex by observation (e.g., trial, trialtype) data matrices of coefficients from first-level GLM.
##    one matrix per subj*wave*session*(run*task)*roi.
##    saved in HDF5 format (.h5) within hidden directory out/parcellated/.d
##  - data_info file in out/parcellated that contains file_name and dset_name for each dataset
## Notes:
##  - removes n-dim B-coefficient vectors (where n is n_vertex) that are all zero (corresponding to all-zero columns 
##    in the time-series design matrix)
## --------------------------------------------------------------------------------------------------------------------


## setup ----

variables <- c("target_denoise", "distractor_denoise", "congruency_denoise")
fnames <- c(hrf = "betas_fithrf.h5", hrf_rr = "betas_rr.h5")# hrf_denoise, hrf_denoise_rr
glmname <- "glmsingle_wave1"
atlas_name <- "glasser2016"
roi_col <- "parcel"
space <- "fsaverage5"
subjlist <- "wave1_unrel_pilot"
sessions <- c("baseline", "proactive", "reactive")
n_cores <- 5
overwrite <- TRUE
delete_files <- FALSE

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
library(mfutils)
source(here("src", "stroop-rsa-pc.R"))

## set variables

task <- "Stroop"
subjects <- as.character(fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]])
tinfo <- read_trialinfo()[subj %in% subjects & wave %in% waves & session %in% sessions]
setkey(tinfo, subj, wave, session, run)  ## for quick subsetting within loop below
atlas <- load_atlas(atlas_name, space)
roiset <- paste0(atlas_name, "_", roi_col)
rois <- unique(atlas$key[[roi_col]])
waves <- c("wave1")

## execute ----

cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
res <- 
    foreach(subj_i = seq_along(subjects), .inorder = FALSE, .combine = "c") %:% 
    foreach(wave_i = seq_along(waves), .inorder = FALSE, .combine = "c") %dopar% {

    for (variable_i in seq_along(variables)) {  ## run these in serial b/c going in same hdf5 file
    for (fname_i in seq_along(fnames)) {
        
        sub <- subjects[subj_i]
        wav <- waves[wave_i]
        variable <- variables[variable_i]
        fname <- fnames[fname_i]
        
        file_name <- here("out", "glms", sub, wav, "RESULTS", "Stroop", glmname, variable, fname)
        for (ses in sessions) {

            tlabels <- tinfo[.(sub, wav, ses), item]  ## get trialtypes
            parcs <- 
                h5read(file_name, ses) %>%
                rename_dim(tlabels) %>% ## add trialtypes to colnames
                .[, Var(., 2, na.rm = TRUE) != 0L] %>% ## filter regressors that are all zero
                parcellate(atlas, roi_col)  ## split matrix into list of length n_roi

            for (i in seq_along(parcs)) {
                for (run_i in 1:2){
                    
                    run <- paste0("run", run_i)
                    run_idx <- 1:n_trial[ses] + n_trial[ses]*(run_i-1)
                    parcrun_i <- parcs[[i]][, run_idx]

                    write_dset(
                        mat = parcrun_i, 
                        roi = names(parcs)[i], 
                        dset_prefix = "coefs",
                        subject = sub, wave = wav, session = ses, run = run, 
                        roiset = roiset, glmname = paste0(glmname, "_", variable, "_", names(fname)), prewh = "none",
                        write_colnames = TRUE, delete_file = delete_files
                        )
                }
            }
        }
        
    }
    }
}

stopCluster(cl)
