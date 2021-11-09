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
    outfile_prefix <- "trr-meanpatt__subjlist-ispc_retest"
    subjects <- fread(here("out/subjlist_all_retest.txt"))[[1]]
    waves <- c("wave1", "wave2")
    sessions <- c("reactive", "proactive", "baseline")
    ii <- 1
    wave_i = 1
} else {
    args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
    glmname <- args$glmname
    roiset <- args$roiset
    prewh <- args$prewh
    outfile_prefix <- args$outfile_prefix
    subjects <- fread(here(args$subjlist))[[1]]
    waves <- strsplit(args$waves, " ")[[1]]
    sessions <- strsplit(args$sessions, " ")[[1]]
}
atlas <- read_atlas(roiset)
rois <- names(atlas$roi)
filename_master <- here("out", "parcellated", "master.h5")
info_master <- read_master(filename_master)
info_master <- info_master[
    prefix == "coefs" & glm %in% paste0("glm-", glmname) & roiset == roiset & prewh %in% prewh &
    subj %in% subjects & wave %in% waves & session %in% sessions
    ]

ttypes_congr <- c(paste0(colors_bias, toupper(colors_bias)), paste0(colors_pc50, toupper(colors_pc50)))

## execute ----

info_master[, g := paste0(subj, "__", session, "__", roi)]
l <- split(info_master, by = "g")  ## NB: this calls split.data.table NOT split.data.frame! must use "by" arg.

cl <- makeCluster(n_core - 2, type = "FORK")
registerDoParallel(cl)
res <- foreach(ii = seq_along(l), .final = function(x) setNames(x, names(l))) %dopar% {
    
    input_val <- l[[ii]]
    sub <- unique(input_val$subj)
    ses <- unique(input_val$session)
    roi <- unique(input_val$roi)

    stopifnot(nrow(input_val) == 4)  ## one for each run*wave

    B <- enlist(waves)
    for (wave_i in seq_along(waves)) {
        hilo_wave <- enlist(runs)
        for (run_i in seq_along(runs)) {
            b <- read_dset(filename_master, input_val[wave == waves[wave_i] & run == runs[run_i], name])
            b_bar <- average(b, colnames(b))
            congruency <- ifelse(colnames(b_bar) %in% ttypes_congr, "congr", "incon")
            hilo_wave[[run_i]] <- average(b_bar, congruency) %*% rbind(-1, 1)
        }
        B[[wave_i]] <- c(Reduce("+", hilo_wave))
        # B[[wave_i]] <- 
        #     colSums(read_dset(filename_master, input_val[wave == waves[wave_i] & run == "run1", name])) +
        #     colSums(read_dset(filename_master, input_val[wave == waves[wave_i] & run == "run2", name]))
    }
    
    Reduce(cor, B)

}
stopImplicitCluster()


d <- as.matrix(abind(res))
d <- as.data.table(d, keep.rownames = TRUE)
d <- separate(d, rn, c("subj", "session", "roi"), sep = "__")

fout <- paste0(outfile_prefix, "__glm-", glmname, "__roiset-", roiset, "__prewh-", prewh, ".RDS")
saveRDS(d, here("out", "res", fout))
