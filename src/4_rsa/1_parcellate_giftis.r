#!/usr/bin/env bash Rscript --vanilla

## --------------------------------------------------------------------------------------------------------------------
## Script: 
## Author: Mike Freund (mfreundc@gmail.com)
## Input: 
##  - 
##  ...
## Output:
##  - 
##  ...
## Notes:
##  ...
##  ...
##  ...
##  
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
#h5disableFileLocking()

## set variables

prewh <- "none"
task <- "Stroop"

if (interactive()) { 
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018_control"
    subjects <- c("115825", "130518", "155938", "178647")
    sessions <- "reactive"
    waves <- c("wave1", "wave2")
    glm_i <- 1
} else {
    args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
    glmname <- args$glmname
    roiset <- args$roiset
    subjects <- fread(here(args$subjlist))[[1]]
    waves <- strsplit(args$waves, " ")[[1]]
    sessions <- strsplit(args$sessions, " ")[[1]]
    if (length(glmname) != 1L || length(roiset) != 1L) stop("not configured for multiple glms or roisets")
}


## read triallabels and atlas

#b <- fread(here("in", "behavior-and-events_stroop_2021-10-20_nice.csv"))
#b <- b[subj %in% subjects & session %in% sessions & wave %in% waves]
#setkey(b, subj, wave, session, run, trial_num)  ## sort to match col order of fmri beta matrix
atlas <- read_atlas(roiset)
rois <- names(atlas$roi)

## get gifti filenames (builds data.table of input filenames, one per hemi)

input <- construct_filenames_gifti(subject = subjects, wave = waves, session = sessions, run = runs, glmname = glmname)

## create hdf5 file and group tree

filename_h5 <- construct_filename_h5(prefix = "coef", glmname = glmname, roiset = roiset, prewh = prewh)
if (!file.exists(filename_h5)) {
    was_created <- h5createFile(filename_h5)
    create_h5groups(filename_h5, subjects = subjects, waves = waves, sessions = sessions, runs = runs, rois = rois)
}


## execute ----


l <- split(input, interaction(input$subject, input$wave, input$session, input$run))

cl <- makeCluster(n_core - 2, type = "FORK")
registerDoParallel(cl)
foreach(glm_i = seq_along(l)) %dopar% {
    
    input_val <- l[[glm_i]]
    subject <- unique(input_val$subj)
    wave <- unique(input_val$wave)
    session <- unique(input_val$session)
    run <- unique(input_val$run)
    
    giftis <- lapply(input_val$filename, read_gifti)
    parcs <- 
        giftis %>%
        concat_hemis(input_val$hemi, pattern = "_Coef") %>% ## extract data and concatenate
        parcellate_data(atlas)  ## split matrix into list of length n_roi
    #triallabs <- b[subj == subject & wave == wave & session == session & run == run]$item  ## get trialtypes
    data_names <- 
        get_path_h5(subject = subject, wave = wave, session = session, run = run, roi = rois, dataname = "coefs")
    ## write:
    Map(function(x, y, file) h5write(obj = x, name = y, file = file), parcs, data_names, file = filename_h5)
    h5closeAll()
    #h5write(parcs[[1]], data_names[[1]], file = filename_h5)
}
stopImplicitCluster()

file_summary <- as.data.table(h5ls(filename_h5))
fwrite(file_summary, gsub(".h5$", ".csv", filename_h5))