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

## glm dir and subject list, misc strings

dir_analysis <- here("out", "glms")
s <- fread(here("out", "subjlist.csv"))
source(here("src", "stroop-rsa-pc.R"))
source(here("src", "_constants.R"))

## command args

if (interactive()) { 

    subjs <- unique(s[is_ispc_retest == TRUE]$subj)[1:2]
    glms <- "lsall_1rpm"
    tasks <- "Stroop"
    sessions <- "reactive"
    waves <- "wave1"
    file_prefix <- "STATS"  ## must be length one
    file_suffix <- "REML.func.gii"
    roiset <- "schaefer"
    result_i <- 1
    run_i <- 1
    roi_i <- 1

} else {

    args <- commandArgs(trailingOnly = TRUE)

}

atlas <- read_atlas(roiset)
names_core32 <- atlas$key$parcel[core32]
rois <- split(atlas$key$parcel, atlas$key$parcel)[names_core32]
rois <- c(rois, Vis = list(atlas$key$parcel[grep("_Vis_", atlas$key$parcel)]), SomMot = list(atlas$key$parcel[grep("_SomMot_", atlas$key$parcel)]))



## build data.table of input arguments

input <- as.data.table(expand.grid(subj = subjs, task = tasks, wave = waves, session = sessions, glm = glms, run = runs))



## execute ----


dir_results <- file.path(dir_analysis, input$subj, input$wave, "RESULTS", input$task, paste0(input$session, "_", input$glm))

for (result_i in seq_along(dir_results)) {
    
    ## read giftis

    dir_result <- dir_results[result_i]

    files <- file.path(dir_result, paste0(file_prefix, "_", input[result_i]$subj, "_", run_hemis, "_", file_suffix))  ## L then R
    gii <- lapply(files, read_gifti)
    parcim <- parcellated_image(gii, fold_hemis = run_hemis, atlas = atlas, pattern = "_Coef")  ## construct parcellated image

    ## save in hdf5
    
    
}

