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

if (interactive()) {  ## add variables (potentially unique to this script) useful for dev
    glmname <- "lsall_1rpm"
    atlas_name <- "glasser2016"
    roi_col <- "superparcel"
    space <- "fsaverage5"
    subjlist <- "baseline_wave1"
    subjects <- fread(here("out", "subjs", paste0("subjlist_", subjlist, ".txt")))[[1]]
    sessions <- "baseline"#"reactive"#c("baseline", "proactive", "reactive")
    waves <- c("wave1")
    glm_i <- 94
    n_cores <- 20
    overwrite <- FALSE
    delete_files <- FALSE
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

if (glmname == "lssep_1rpm") {
    suffix <- ".gii"
    pattern <- NULL
} else if (glmname %in% c("lsall_1rpm", "condition_1rpm", "lsall_acompcor06_1rpm")) {
    suffix <- "_REML.func.gii"
    pattern <- "_Coef"
} else {
    stop("not configured for provided glmname")
}
tinfo <- read_trialinfo()[subj %in% subjects & wave %in% waves & session %in% sessions]
setkey(tinfo, subj, wave, session, run)  ## for quick subsetting within loop below
atlas <- load_atlas(atlas_name, space)
atlas <- add_superparcels(atlas, superparcels, roi_col = roi_col)
roiset <- paste0(atlas_name, "_", roi_col)
rois <- unique(atlas$key[[roi_col]])
rois <- rois[!is.na(rois)]

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
    if (glmname == "condition_1rpm") tlabels <- c(ttypes$bias, ttypes$pc50)
        
    giftis <- lapply(input_val$filename, read_gifti)
    parcs <- 
        giftis %>%
        concat_hemis(input_val$hemi, pattern = pattern) %>% ## extract data and concatenate
        rename_dim(tlabels) %>% ## add trialtypes to colnames
        .[, Var(., 2) != 0L] %>%  ## filter regressors that are all zero
        parcellate(atlas, roi_col) %>% ## split matrix into list of length n_roi
        .[lengths(.) != 0]  ## remove "parcels" with roi_col set to NA

    out_names <- construct_filename_h5(
        base_dir = here::here("out", "parcellated", ".d2"),
        subject = sub, wave = wav, session = ses, run = run, roiset = roiset, roi = rois
        )
    if (!overwrite) {
        is_extant <- file.exists(out_names)
        if (all(is_extant)) return(list("All parcellated data files already exist. Skipping..."))
        parcs <- parcs[!is_extant]
    }

    nms <- enlist(names(parcs))
    for (i in seq_along(parcs)) {
        nms[[i]] <- write_dset(
            mat = parcs[[i]], roi = names(parcs)[i],
            base_dir = here::here("out", "parcellated", ".d2"), 
            dset_prefix = "coefs",
            subject = sub, wave = wav, session = ses, run = run, 
            roiset = roiset, glmname = glmname, prewh = "none",
            write_colnames = TRUE, 
            delete_file = delete_files
            )
        #print(i)
    }

    nms

}
stopCluster(cl)

data_info <- as.data.table(do.call(rbind, res))
fn <- construct_filename_datainfo(
        prefix = "coefs", subjlist = subjlist, glmname = glmname, roiset = roiset, prewh = "none"
        )
fwrite(data_info, fn)
