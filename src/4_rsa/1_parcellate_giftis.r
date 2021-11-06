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


## execute ----


l <- split(input, interaction(input$subject, input$wave, input$session, input$run))

cl <- makeCluster(n_core - 2, type = "FORK")
registerDoParallel(cl)
names_for_link <- foreach(glm_i = seq_along(l), .inorder = FALSE, .combine = "c") %dopar% {
    
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
    nms <- Map(
        function(x, y, ...) write_dset(mat = x, roi = y, ...), 
        parcs, names(parcs),
        dset_prefix = "coefs",
        subject = subject, wave = wave, session = session, run = run, 
        roiset = roiset, glmname = glmname, prewh = prewh
        )
    nms
}
stopImplicitCluster()

## update master.h5
write_links(names_for_link)
#fid <- H5Fopen(here("out", "parcellated", "master.h5"))
#h5closeAll()
