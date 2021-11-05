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


## todo
## -what to do about subjects input? make separate file for each subject group and input as file?
## -roiset to functions
## -tests for parcellated_image()


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
source(here("src", "stroop-rsa-pc.R"))

## set variables

prewh <- "none"
task <- "Stroop"

if (interactive()) { 
    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018_control"
    subjects <- c("115825", "132017")
    sessions <- "reactive"
    waves <- c("wave1", "wave2")
} else {
    args <- commandArgs(trailingOnly = TRUE)
}

## read triallabels and atlas

b <- fread(here("in", "behavior-and-events_stroop_2021-10-20_nice.csv"))
b <- b[subj %in% subjects & session %in% sessions & wave %in% waves]
setkey(b, subj, wave, session, run, trial_num)  ## sort to match col order of fmri beta matrix
atlas <- read_atlas(roiset)
rois <- names(atlas$roi)

## get gifti filenames (builds data.table of input filenames, one per hemi)

input <- construct_filenames_gifti(subject = subjects, wave = waves, session = sessions, run = runs, glmname = glmname)

## create hdf5 file and group tree

filename_h5 <- construct_filename_h5(prefix = "coef", glmname = glmname, roiset = roiset, prewh = prewh)
was_created <- h5createFile(filename_h5)
create_h5groups(filename_h5, subjects = subjects, waves = waves, sessions = sessions, runs = runs, rois = rois)



## execute ----

l <- split(input, interaction(input$subject, input$wave, input$session, input$run))

for (glm_i in seq_along(l)) {
    
    input_val <- l[[glm_i]]

    ## extract info
    
    subject <- unique(input_val$subj)
    wave <- unique(input_val$wave)
    session <- unique(input_val$session)
    run <- unique(input_val$run)
    
    giftis <- lapply(input_val$filename, read_gifti)
    parcs <- 
        giftis %>%
        concat_hemis(input_val$hemi, pattern = "_Coef") %>% ## extract data and concatenate
        parcellate_data(atlas)  ## split matrix into list of length n_roi
    triallabs <- b[subj == subject & wave == wave & session == session & run == run]$item  ## get trialtypes
    ## store trialtypes in parcs?
    ## save in hdf5

}

