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

## read behavioral data and subject list
b <- fread(here("in", "behavior-and-events_stroop_2021-10-20_nice.csv"))
s <- fread(here("out", "subjlist.csv"))


## set variables

prewh <- "none"
task <- "Stroop"

if (interactive()) { 

    glmname <- "lsall_1rpm"
    roiset <- "Schaefer2018_control"
    
    subjects <- unique(s[is_ispc_retest == TRUE]$subj)[1:2]
    sessions <- "reactive"
    waves <- c("wave1", "wave2")
    file_prefix <- "STATS"  ## must be length one
    file_suffix <- "REML.func.gii"
    do_dev_rois <- TRUE
    result_i <- 1
    run_i <- 1
    roi_i <- 1

} else {

    args <- commandArgs(trailingOnly = TRUE)

}

## read atlas

atlas <- read_atlas(roiset)

## wrangle behavioral data
b <- b[subj %in% subjects & session %in% sessions & wave %in% waves]
b <- arrange(b, subj, wave, session, run, trial_num)

## build data.table of input arguments
## each row corresponds to a single subj*wave*session*run directory, which contains two giftis (L, R hemi)
input <- as.data.table(expand.grid(
    subj = subjs, wave = waves, session = sessions, run = runs,
    glm = glm_name,
    stringsAsFactors = FALSE
    )
)
input$dir_results <- file.path(
    dir_analysis, input$subj, input$wave, "RESULTS", task, paste0(input$session, "_", input$glm)
    )

## execute ----

fname <- here(
    "out", 
    construct_filename(
        data_type = "pimage", prefix = "B", filetype = ".h5", 
        glm_name = glm_name, roi_set = roiset, prewhitened = prewhitened
        )
)

was_created <- h5createFile(fname)

for (result_i in seq_len(nrow(input))) {
    
    ## extract info
    dir_result <- input[result_i]$dir_result
    subject <- input[result_i]$subj
    session <- input[result_i]$session
    wave <- input[result_i]$wave
    
    ## read giftis (L then R)
    files <- file.path(dir_result, paste0(file_prefix, "_", input[result_i]$subj, "_", run_hemis, "_", file_suffix))
    giftis <- lapply(files, read_gifti)
    pimage <- giftis_to_parcellated_image(
        giftis, 
        folds = gsub("_L|_R", "", run_hemis), 
        hemis = gsub("run[0-9]_", "", run_hemis),
        atlas = atlas, 
        colname_data = "B",
        subject = subject, session = session, wave = wave, task = task,
        glm_name = glm_name, roi_set = roiset, prewhitened = "none", 
        shrinkage_var = as.numeric(NA), shrinkage_cov = as.numeric(NA),
        pattern = "_Coef"
        )  ## construct parcellated image
    
    ## add trialtype info
    b_val <- b[subj == input[result_i]$subj & wave == input[result_i]$wave & session == input[result_i]$session]  ## includes both runs
    trialtypes <- list(run1 = b_val[run == 1]$item, run2 = b_val[run == 2]$item)
    #pimage <- relabel(pimage, "data", trialtypes)
    
    ## save in hdf5

    create_nested_group(g = c(subject, session, wave, unique(pimage$data$roi), fold), file_name = fname)
    write_nested_h5(pimage, colname_data = "B", file_name = fname)

}





# colname_data = "B"
# file_name = fname_h5
# data_name = colname_data
# file_i <- 1


## TODO:
## writing a file into database
## - subj/session/wave/roi/run: data, bad_verts, labels, meta -- in filename
## reading a file from database and building a parcellated_data/image object
## - specify filename: subj/session/wave/roi/run
## - generate metadata and data.table structure
## - read in matrices and set to data.table
## constructor and validation

## (lsa, lss) * (parcel, network, rois) * (prew_resid, prew_obs, prew_rest, none)

## largest array possible:
## - session, subj, wave, roi: feature, obs, fold
## what may be added?

#metadata to group/dims
## labels; rm row/colnames
## add rownames function
