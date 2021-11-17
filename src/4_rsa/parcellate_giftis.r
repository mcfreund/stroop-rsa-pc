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

task <- "Stroop"

if (interactive()) {  ## add variables (potentially unique to this script) useful for dev
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018Dev"
    subjlist <- "ispc_retest"
    subjects <- fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]][1:5]
    sessions <- "reactive"
    waves <- c("wave1", "wave2")
    glm_i <- 1
    n_cores <- 10
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

if (glmname == "lssep_1rpm") {
    suffix <- ".gii"
    pattern <- NULL
} else if (glmname == "lsall_1rpm") {
    suffix <- "_REML.func.gii"
    pattern <- "_Coef"
}
tinfo <- read_trialinfo()[subj %in% subjects & wave %in% waves & session %in% sessions]
setkey(tinfo, subj, wave, session, run)  ## for quick subsetting within loop below
atlas <- read_atlas(roiset)
rois <- names(atlas$roi)


## execute ----

input <- construct_filenames_gifti(
    subject = subjects, wave = waves, session = sessions, run = runs, glmname = glmname,
    suffix = suffix
    )
input[, g := paste0(subject, "__", wave, "__", session, "__", run)]
l <- split(input, by = "g")

cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
res <- foreach(glm_i = seq_along(l), .inorder = FALSE, .combine = "c") %dopar% {
#for (glm_i in seq_along(l)) {
    
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
        concat_hemis(input_val$hemi, pattern = pattern) %>% ## extract data and concatenate
        rename_dim(tlabels) %>%  ## add trialtypes to colnames
        parcellate_data(atlas) ## split matrix into list of length n_roi
    nms <- Map(
        function(x, y, ...) write_dset(mat = x, roi = y, ...), 
        parcs, names(parcs),
        dset_prefix = "coefs",
        subject = sub, wave = wav, session = ses, run = run, 
        roiset = roiset, glmname = glmname, prewh = "none",
        write_colnames = TRUE
        )

    nms

}
stopImplicitCluster()

data_info <- as.data.table(do.call(rbind, res))
fn <- construct_filename_datainfo(
        prefix = "coefs", subjlist = subjlist, glmname = glmname, roiset = roiset, prewh = "none"
        )
fwrite(data_info, fn)
