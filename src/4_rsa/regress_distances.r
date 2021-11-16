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

if (interactive()) {  ## add variables (potentially unique to this script) useful for dev
    glmname <- "ls_1rpm"
    roiset <- "Schaefer2018Dev"
    prewh <- "none"
    measure <- "cveuc"  ## "crcor"
    subjlist <- "ispc_retest"
    subjects <- fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]][1:5]
    waves <- c("wave1", "wave2")
    sessions <- "reactive"
    ttype_subset <- "bias"
    subject = subjects[1]
    wave = waves[1]
    subject_i = 1
    wave_i = 1
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

stopifnot(sessions %in% c("baseline", "proactive", "reactive"))
stopifnot(measure %in% c("crcor", "cveuc"))

atlas <- read_atlas(roiset)
rois <- names(atlas$roi)
X <- read_model_xmat(measure, sessions, ttype_subset)


## execute ----
## loop over sess, measure, prew, ...


D <- read_rdms(
    .measure = measure, .glmname = glmname, .roiset = roiset, .prewh = prewh, 
    .subjects = subjects, .session = sessions, .waves = waves,
    .ttypes1 = intersect(ttypes_by_run[[sessions]]$run1, ttypes[[ttype_subset]]),
    .ttypes2 = intersect(ttypes_by_run[[sessions]]$run2, ttypes[[ttype_subset]])
    )


res <- enlist(combo_paste(subjects, waves, sep = "_"))
for (subject_i in seq_along(subjects)) {
    for (wave_i in seq_along(waves)) {
        
        Di <- D[, , , subject_i, wave_i]
        
        ## preproc:
        if (measure == "cveuc") {
            d <- apply(Di, 3, vec)  ## vectorize
            Y <- d %*% diag(1/sqrt(Var(d, 2)))  ## scale
        } else if (measure == "crcor") {
            d <- apply(Di, 3, c)  ## vectorize
            Y <- atanh(d)  ## Fisher's z transform
        }
        
        B <- coef(.lm.fit(x = X, y = Y))
        
        ## save:
        nm <- paste0(subjects[subject_i], "_", waves[wave_i])
        res[[nm]] <- B %>% 
            tidy_model(terms = colnames(X), outcomes = dimnames(Di)$roi) %>% 
            melt(id.vars = "term", value.name = "b", variable.name = "roi")
        
    }
}

data <- rbindlist(res, idcol = "subject_wave")
data <- separate(data, subject_wave, into = c("subject", "wave"))

fout <- paste0("weights-", measure, "__subjlist-", subjlist, "__glm-", glmname, "__roiset-", roiset, "__prewh-", prewh, ".csv")
fwrite(data, here("out", "res", fout))
