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

# l <- split(input, interaction(input$subject, input$wave, input$session, input$run))

# cl <- makeCluster(n_core - 2, type = "FORK")
# registerDoParallel(cl)
# names_for_link <- foreach(glm_i = seq_along(l), .inorder = FALSE, .combine = "c") %dopar% {
    
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
