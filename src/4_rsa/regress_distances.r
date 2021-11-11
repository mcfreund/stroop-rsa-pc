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
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018Dev"
    prewh <- "none"
    measure <- "cveuc"  ## "crcor"
    subjects <- fread(here("out/subjlist_ispc_retest.txt"))[[1]][1:5]
    waves <- c("wave1", "wave2")
    session <- "reactive"
    subject = subjects[1]
    wave = waves[1]
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

atlas <- read_atlas(roiset)
rois <- names(atlas$roi)

tidy_model <- function(B, models, outcomes) {
    B <- as.data.table(B)
    names(B) <- outcomes
    B$model <- models
    B
 }


## execute ----
## loop over sess, measure, prew, ...

if (measure == "cveuc") {
    X <- cbind(x0 = 1, read_model_rdm(cells = "lowertri", session = session))
} else if (measure == "crcor") {
    stop("not configured yet")
}

D <- read_rdms(
    .measure = measure, .glmname = glmname, .roiset = roiset, .prewh = "none", 
    .subjects = subjects, .session = session, .waves = waves
    )
D1 <- read_rdms(
    .measure = measure, .glmname = glmname, .roiset = roiset, .prewh = "obsall", 
    .subjects = subjects, .session = session, .waves = waves
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
            d <- apply(Di, 3, offdiag)  ## vectorize
            Y <- atanh(d)  ## Fisher's z transform
        }
        
        B <- coef(.lm.fit(x = X, y = Y))
        
        ## save:
        nm <- paste0(subjects[subject_i], "_", waves[wave_i])
        res[[nm]] <- B %>% 
            tidy_model(models = colnames(X), outcomes = dimnames(Di)$roi) %>% 
            melt(id.vars = "model", value.name = "b", variable.name = "roi")
        
    }
}

data <- rbindlist(res, idcol = "subject_wave")
data <- separate(data, subject_wave, into = c("subject", "wave"))

fout <- paste0("weights__glm-", glmname, "__roiset-", roiset, "__prewh-", prewh, ".csv")
fwrite(data, here("out", "res", fout))
