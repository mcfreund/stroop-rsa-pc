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
##  - master.h5 file in out/parcellated that contains external (soft) links to each h5 dataset in .d/
## Notes:
##  - dataset links in master.h5 file cannot be modified once created?
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
    roiset <- "Schaefer2018Dev"
    subjects <- c("115825", "130518", "155938", "178647")
    #subjects <- fread(here("out/subjlist_ispc_retest.txt"))[[1]]
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

tinfo <- read_trialinfo()[subj %in% subjects & wave %in% waves & session %in% sessions]
setkey(tinfo, subj, wave, session, run)  ## for quick subsetting within loop below
atlas <- read_atlas(roiset)
rois <- names(atlas$roi)
input <- construct_filenames_gifti(subject = subjects, wave = waves, session = sessions, run = runs, glmname = glmname)


## execute ----

input[, g := paste0(subject, "__", wave, "__", session, "__", run)]
l <- split(input, by = "g")

cl <- makeCluster(n_core - 2, type = "FORK")
registerDoParallel(cl)
names_for_link <- foreach(glm_i = seq_along(l), .inorder = FALSE, .combine = "c") %dopar% {
    
    input_val <- l[[glm_i]]
    sub <- unique(input_val$subj)
    wav <- unique(input_val$wave)
    ses <- unique(input_val$session)
    run <- unique(input_val$run)
    run_i <- as.numeric(gsub("run", "", run))
    tlabels <- tinfo[.(sub, wav, ses, run_i), item]  ## get trialtypes
    
    giftis <- lapply(input_val$filename, read_gifti)
    parcs <- 
        giftis %>%
        concat_hemis(input_val$hemi, pattern = "_Coef") %>% ## extract data and concatenate
        rename_dim(tlabels) %>%  ## add trialtypes to colnames
        parcellate_data(atlas) ## split matrix into list of length n_roi
    nms <- Map(
        function(x, y, ...) write_dset(mat = x, roi = y, ...), 
        parcs, names(parcs),
        dset_prefix = "coefs",
        subject = sub, wave = wav, session = ses, run = run, 
        roiset = roiset, glmname = glmname, prewh = prewh,
        write_colnames = TRUE
        )

    nms

}
stopImplicitCluster()


write_links(names_for_link)  ## update master.h5
