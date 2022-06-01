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
library(mfutils)
source(here("src", "stroop-rsa-pc.R"))


## set variables

task <- "Stroop"

if (interactive()) {  ## add variables (potentially unique to this script) useful for dev
    glmname <- "lsall_1rpm"
    prewh <- "none"
    atlas_name <- "schaefer2018_17_200"
    roi_col <- "parcel"
    space <- "fsaverage5"
    measure <- "crcor"
    subjlist <- "proactive_wave2"
    subjects <- fread(here("out", "subjs", paste0("subjlist_", subjlist, ".txt")))[[1]]
    waves <- "wave2"
    sessions <- "proactive"
    #sessions <- "baseline"
    ttype_subset <- "bias"
    subject = subjects[1]
    suffix <- "__seswave-baseline_wave1"
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

stopifnot(sessions %in% c("baseline", "proactive", "reactive"))
stopifnot(measure %in% c("crcor", "cveuc"))
if (!exists("suffix")) suffix <- ""
subjects <- as.character(subjects)

atlas <- load_atlas(atlas_name, space)
atlas <- add_superparcels(atlas, superparcels, roi_col = roi_col)
roiset <- paste0(atlas_name, "_", roi_col)
rois <- unique(atlas$key[[roi_col]])
rois <- rois[!is.na(rois)]
## remove hippocampi from glasser atlas (not represented in fsaverages???)
if (atlas_name == "glasser2016" && grepl("fsaverage", space)) {
    rois <- rois[rois != "L_H" & rois != "R_H"]
}
roiset <- paste0(atlas_name, "_", roi_col)

X <- enlist(sessions)
for (session in sessions) X[[session]] <- read_model_xmat(measure, session, ttype_subset)


## execute ----
## loop over sess, measure, prew, ...

D <- enlist(combo_paste(sessions, "_", waves))
for (session in sessions) {
    for (wave in waves) {
        D[[paste0(session, "_", wave)]] <- 
            read_rdms(
                .measure = measure, .prewh = paste0(prewh, "__2022-05-31"),  ## need to change this manually....
                .glmname = glmname, .roiset = roiset, .waves = wave,
                .session = session, .subjects = subjects,
                .ttype_subset = ttype_subset,
                .ttypes1 = intersect(ttypes_by_run[[session]]$run1, ttypes[[ttype_subset]]),
                .ttypes2 = intersect(ttypes_by_run[[session]]$run2, ttypes[[ttype_subset]]),
                .rois = rois
            )
    }
}

res <- enlist(combo_paste(subjects, waves, sessions, sep = "_"))
for (subject_i in seq_along(subjects)) {
    for (wave_i in seq_along(waves)) {
        for (session_i in seq_along(sessions)) {

        session <- sessions[session_i]
        wave <- waves[wave_i]
        subject <- subjects[subject_i]

        Di <- D[[paste0(session, "_", wave)]][, , , subject_i, wave_i]

        ## preproc:
        if (measure == "cveuc") {
            d <- apply(Di, 3, squareform)  ## vectorize
            Y <- d %*% diag(1/sqrt(Var(d, 2)))  ## scale
        } else if (measure == "crcor") {
            d <- apply(Di, 3, squareform, lt = FALSE)  ## vectorize
            Y <- atanh(d)  ## Fisher's z transform
        }
        
        B <- coef(.lm.fit(x = X[[session]], y = Y))
        
        ## save:
        nm <- paste0(subject, "_", wave, "_", session)
        res[[nm]] <- B %>% 
            tidy_model(terms = colnames(X[[session]]), outcomes = dimnames(Di)$roi) %>% 
            melt(id.vars = "term", value.name = "b", variable.name = "roi")
        
        }
    }
}

data <- rbindlist(res, idcol = "subject_wave_session")
data <- separate(data, subject_wave_session, into = c("subject", "wave", "session"))

fout <- construct_filename_weights(
    measure = measure, subjlist = subjlist, glmname = glmname, roiset = roiset, ttype_subset = ttype_subset, 
    prewh = prewh, suffix = suffix,
    base_dir = here("out", "weights")
    )
fwrite(data, fout)
