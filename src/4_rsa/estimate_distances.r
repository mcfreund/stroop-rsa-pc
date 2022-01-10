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
    space <- "fsaverage5"
    roi_col <- "parcel"
    subjlist <- "mi1"
    subjects <- fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]]
    waves <- "wave1"
    sessions <- "proactive"
    measure <- "crcor"  ## crcor
    prewh <- "none"
    ttype_subset <- "all"
    ii <- 596
    n_cores <- 28
    run_i <- 1
    n_resamples <- 1E4
    overwrite <- TRUE
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

stopifnot(sessions %in% c("baseline", "proactive", "reactive"))
stopifnot(measure %in% c("crcor", "cveuc"))

atlas <- load_atlas(atlas_name, space)
rois <- unique(atlas$key[[roi_col]])
## remove hippocampi from glasser atlas (not represented in fsaverages???)
if (atlas_name == "glasser2016" && grepl("fsaverage", space)) {
    rois <- rois[rois != "L_H" & rois != "R_H"]
}
roiset <- paste0(atlas_name, "_", roi_col)


## execute ----


input <- construct_filenames_h5(
    prefix = "coefs", subjects = subjects, waves = waves, sessions = sessions, rois = rois, runs = runs, 
    glmname = glmname, prewh = prewh
)
input[, g := paste0(subj, "__", wave, "__", session, "__", roi)]
l <- split(input, by = "g")


## skip computation of existing output files if overwrite == FALSE:

fname <- construct_filename_rdm(
    measure = measure, glmname = glmname, ttype_subset = ttype_subset, roiset = roiset, prewh = prewh
    )  ## output file

if (!overwrite) {
    if (file.exists(fname)) {
        group_name <- combo_paste(subjects, waves, sessions, sep = "__")
        group_exists <- group_name %in% h5ls(fname)$name
        if (any(group_exists)) {
            existing_group_idx <- grep(paste0(group_name[group_exists], collapse = "|"), names(l))
            l <- l[-existing_group_idx]
            if (length(l) == 0L) stop("All to-be-written files already exist. Exiting...")
        }
    }
}


## precompute index matrices for crcor resampling:

if (measure == "crcor") {
    pick_an_roi <- unique(input$roi)[1]  ## given subj*sess*wave*run, roi doesn't matter b/c cond orders all same
    is_roi <- vapply(l, function(x) unique(x$roi) == pick_an_roi, logical(1))
    l_idx <- l[is_roi]
    names(l_idx) <- gsub(paste0("__", pick_an_roi), "", names(l_idx))  ## remove roi name from list names
    l_idx <- mclapply(
        l_idx, 
        function(x) {
            ses <- unique(x$session)
            nms1 <- colnames(read_dset(x[run == "run1", file_name], x[run == "run1", dset_name]))  ## cond orders run 1
            nms2 <- colnames(read_dset(x[run == "run2", file_name], x[run == "run2", dset_name]))  ## "  " run2
            ## pull only trialtypes being analyzed:
            nms1 <- nms1[nms1 %in% intersect(ttypes_by_run[[ses]]$run1, ttypes[[ttype_subset]])]
            nms2 <- nms2[nms2 %in% intersect(ttypes_by_run[[ses]]$run2, ttypes[[ttype_subset]])]
            ## resample trial indices:
            idx1 <- get_resampled_idx(nms1, n_resamples, seed = 0)
            idx2 <- get_resampled_idx(nms2, n_resamples, seed = 0)
            list(run1 = idx1, run2 = idx2)
        },
        mc.cores = n_cores
    )
}


cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
res <- foreach(
    ii = seq_along(l[1:4]), 
    .final = function(x) setNames(x, names(l)[1:4]), 
    .noexport = "crcor_rcpp"
    ) %dopar% {
#resall <- enlist(names(l))
#for (ii in seq_along(l)) {
#for (ii in seq(ii, length(l))) {

    input_val <- l[[ii]]
    ses <- unique(input_val$session)
    
    ## read betas and subset trialtypes, vertices
    
    B <- enlist(runs)
    is_bad_vert <- enlist(runs)
    for (run_i in seq_along(runs)) {
        Bi <- read_dset(input_val[run == runs[run_i], file_name], input_val[run == runs[run_i], dset_name])
        ttypes_to_get <- intersect(ttypes_by_run[[ses]][[run_i]], ttypes[[ttype_subset]])
        Bi <- Bi[, colnames(Bi) %in% ttypes_to_get]  ## discard low-trialcount trialtypes
        B[[run_i]] <- Bi
        is_bad_vert[[run_i]] <- is_equal(Var(Bi, 1), 0)  ## verts with no BOLD
    }
    ## subset vertices:
    is_good_vert <- !Reduce("|", is_bad_vert)
    if (sum(is_good_vert) == 0) {
        stop(c("no good verts: ", paste0(input_val[1, 1:5], sep = " ")))
    } else {
        B <- lapply(B, function(x) x[is_good_vert, ])  ## discard vertices with no BOLD
    }
    
    ## estimate distances/similarities

    if (measure == "cveuc") {  ## cross-validated euclidean
        res <- cvdist(
            average(B$run1, g = colnames(B$run1)),
            average(B$run2, g = colnames(B$run2)),
            scale = TRUE
        )
    } else if (measure == "crcor") {    ## cross-run correlation (with downsampling)
        id <- unique(input_val[, paste0(subj, "__", wave, "__", session)])
        res <- crcor(x1 = B$run1, x2 = B$run2, idx1 = l_idx[[id]]$run1, idx2 = l_idx[[id]]$run2)
    }

    #resall[[ii]] <- res
    res

}
stopCluster(cl)


## write arrays within hdf5 file

if (!file.exists(fname)) h5createFile(fname)
outnames <- unique(rbindlist(l)[, .(subj, wave, session)])
for (row_i in seq_len(nrow(outnames))) {

    nm <- paste0(outnames[row_i, ], collapse = "__")

    ## extract data and concatenate into 3D array (RDM, roi)
    dat <- res[grep(nm, names(res))]
    dat <- abind(dat, rev.along = 0)
    dimnames(dat)[[3]] <- gsub(paste0(nm, "__"), "", dimnames(dat)[[3]])
    
    h5write(dat, fname, nm)

}
