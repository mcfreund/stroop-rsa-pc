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
    subjects <- c("115825", "130518", "155938", "178647")
    #subjects <- fread(here("out/subjlist_ispc_retest.txt"))[[1]]
    sessions <- "reactive"
    waves <- c("wave1", "wave2")
    prewh <- "none"
    glm_i <- 1
} else {
    args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
    glmname <- args$glmname
    roiset <- args$roiset
    prewh <- args$prewh
    subjects <- fread(here(args$subjlist))[[1]]
    waves <- strsplit(args$waves, " ")[[1]]
    sessions <- strsplit(args$sessions, " ")[[1]]
    if (length(glmname) != 1L || length(roiset) != 1L || length(prewh) != 1L) {
        stop("not configured for multiple glms, roisets, or prewhs")
    }
}

tinfo <- read_trialinfo()[subj %in% subjects & wave %in% waves & session %in% sessions]
setkey(tinfo, subj, wave, session, run)  ## for quick subsetting within loop below
atlas <- read_atlas(roiset)
rois <- names(atlas$roi)


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

    
#     input_val <- l[[glm_i]]
#     sub <- unique(input_val$subj)
#     wav <- unique(input_val$wave)
#     ses <- unique(input_val$session)
#     run <- unique(input_val$run)
#     run_i <- as.numeric(gsub("run", "", run))
#     tlabels <- tinfo[.(sub, wav, ses, run_i), item]  ## get trialtypes
    
#     giftis <- lapply(input_val$filename, read_gifti)
#     parcs <- 
#         giftis %>%
#         concat_hemis(input_val$hemi, pattern = "_Coef") %>% ## extract data and concatenate
#         rename_dim(tlabels) %>%  ## add trialtypes to colnames
#         parcellate_data(atlas) ## split matrix into list of length n_roi
#     nms <- Map(
#         function(x, y, ...) write_dset(mat = x, roi = y, ...), 
#         parcs, names(parcs),
#         dset_prefix = "coefs",
#         subject = sub, wave = wav, session = ses, run = run, 
#         roiset = roiset, glmname = glmname, prewh = prewh,
#         write_colnames = TRUE
#         )

#     nms

# }
# stopImplicitCluster()


# write_links(names_for_link)  ## update master.h5
