#!/usr/bin/env bash Rscript --vanilla

## input args

if (interactive()) {  ##  for dev:

    glmname <- "lsall_1rpm"
    atlas_name <- "glasser2016"  ## "schaefer2018_17_200"
    space <- "fsaverage5"
    roi_col <- "parcel"
    subjlist <- "wave1_unrel"
    dowave <- 1
    n_cores <- 24
    n_resamples <- 1E2
    decoder <- "svm"
    center <- TRUE
    detrend <- TRUE
    degree <- 1
    demean <- TRUE
    divnorm <- TRUE
    stimset <- "pc50"
    ii = 1
    roi_i <- 1

} else {

    library(argparse)

    parser <- ArgumentParser(description = "Decode target, distractor, and congruency; train bas, test pro/rea.")
    parser$add_argument("atlas_name", type = "character")
    parser$add_argument("roi_col", type = "character")
    parser$add_argument("subjlist", type = "character")
    parser$add_argument("--space", type = "character", default = "fsaverage5")
    parser$add_argument("--n_cores", type = "double", default = 28)
    parser$add_argument("--n_resamples", type = "double", default = 100)
    parser$add_argument("--decoder", type = "character", default = "svm")
    parser$add_argument("--dowave", type = "double", default = 1)
    parser$add_argument("--glmname", type = "character", default = "lsall_1rpm")
    parser$add_argument("--center", action = "store_true")
    parser$add_argument("--detrend", action = "store_true")
    parser$add_argument("--degree", type = "double", default = 1)
    parser$add_argument("--demean", action = "store_false")
    parser$add_argument("--divnorm", action = "store_false")

    args <- parser$parse_args()
    print(args)
    for (i in seq_along(args)) assign(names(args)[i], args[[i]])

}

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
suppressMessages(library(mda))
library(e1071)

source(here("src", "stroop-rsa-pc.R"))


## set variables

subjects <- fread(here("out", paste0("subjlist_", subjlist, ".txt")))[, V1]
waves <- waves[dowave]
task <- "Stroop"
variables <- c("target", "distractor", "congruency")
train_session <- "baseline"
test_sessions <- c("proactive", "reactive")
trialcounts <- fread(here("in", "trialcounts.csv"))
stimsets <- list(
    bias = trialcounts[pc == "bias", unique(stimulus)],
    pc50 = trialcounts[pc == "pc50", unique(stimulus)],
    biaspc50 = trialcounts[, unique(stimulus)]
)

## stimuli to drop from training set:

stimuli_drop <- list(
    run1 = c("blackBLACK", "pinkPINK", "redRED", "whiteWHITE"),  ## items to drop in run 1
    run2 = c("yellowYELLOW", "greenGREEN", "blueBLUE", "purplePURPLE")  ## items to drop in run 2
)
stimuli_to6 <- c("blueBLUE_run1", "purplePURPLE_run1", "redRED_run2", "whiteWHITE_run2")
stimuli_to3 <- c("blueBLUE_run2", "purplePURPLE_run2", "redRED_run1", "whiteWHITE_run1")
stimuli_eachrun <- c("blackBLACK", "greenGREEN", "pinkPINK", "yellowYELLOW")

## atlas info:

atlas <- load_atlas(atlas_name, space)
atlas <- add_superparcels(atlas, superparcels, roi_col = roi_col)
rois <- unique(atlas$key[[roi_col]])
rois <- rois[!is.na(rois)]
## remove hippocampi from glasser atlas (not represented in fsaverages???)
if (atlas_name == "glasser2016" && grepl("fsaverage", space)) {
    rois <- rois[rois != "L_H" & rois != "R_H"]
}
roiset <- paste0(atlas_name, "_", roi_col)


## utilities

summarize_predictions <- function(fit, newdata, decoder) {

    if (decoder == "svm") {
        x <- attr(predict(fit, newdata = newdata, decision.values = TRUE), "decision.values")
    } else {
        stop("not yet configured for anything but svm.")
    }
    
    x

}


## looping info

input <- construct_filenames_h5(
    prefix = "coefs", subjects = subjects, waves = waves, sessions = sessions, rois = rois, runs = runs, 
    glmname = glmname, prewh = "none",
    base_dir = here("out", "parcellated", ".d2")
)
## identify and remove subjects with missing data:
input[, file_exists := file.exists(file_name)]
subjs_missing <- input[file_exists == FALSE, unique(subj)]
print(noquote(paste0(c("removing subjects due to missing data:", subjs_missing), collapse = " ")))
input <- input[!subj %in% subjs_missing]
## make larger chunks for each core by redefining looping variable (have each core loop over ROIs):
## each row of input_subset defines the location of data for a single subject*wave*session*roi (both runs):
input[, g := subj]  ## run subjects in parallel (ROI in serial)
g <- unique(input$g)


## run ----


time_beg <- Sys.time()
cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
out <- foreach(ii = seq_along(g), .inorder = FALSE) %dopar% {

    id <- g[ii]
    inputs <- input[g == id]
    
    res <- enlist(paste0(id, "__", rois))
    for (roi_i in seq_along(rois)) {
        

        ## read betas:
        B <- enlist(combo_paste(sessions, "_", runs))
        for (session_nm in sessions) {
            for (run_nm in runs) {
                sesrun <- paste0(session_nm, "_", run_nm)
                nms <- inputs[roi == rois[roi_i] & run == run_nm & session == session_nm, c("file_name", "dset_name")]
                B[[sesrun]] <- read_dset(nms$file_name, nms$dset_name)
            }
        }
        

        ## subset vertices to common "good" set:
        n_vert <- nrow(B[[1]])
        is_good_vert <- rowSums(is_equal(vapply(B, Var, numeric(n_vert)), 0)) < 1
        B <- lapply(B, function(x) x[is_good_vert, ])  ## discard vertices with no BOLD signal
        

        ## preprocess each session*run
        
        B_prep <- B  ## copy
        if (center) {
            ## regress mean pattern?
            B_prep <- lapply(B_prep, function(betas) {
                betas_i <- betas[, colnames(betas) %in% stimuli_eachrun]
                mu <- average(betas_i, colnames(betas_i))
                mubar <- cbind(scale2unit(rowMeans(mu)))
                fit <- .lm.fit(x = mubar, y = betas)
                resid(fit)
            })
        }
        if (detrend) {
            ## regress time-related nuisance variance?
            B_prep <- lapply(B_prep, function(betas, degree, n_block) {
                trial_num <- seq_len(ncol(betas))
                n <- length(trial_num)/n_block ## trials per block
                blocks <- cbind(
                    block1 = ifelse(trial_num <= n, trial_num, 0),
                    block2 = ifelse(trial_num > n & trial_num <= n*2, trial_num-n, 0),
                    block3 = ifelse(trial_num > n*2, trial_num-n*2, 0)
                )
                x <- apply(blocks, 2, poly, degree = degree)  ## build polynomial regressors
                dim(x) <- c(n*n_block, degree*n_block)  ## reshape
                X <- cbind(1, x)
                fit <- .lm.fit(x = X, y = betas)
                resid(fit)
                },
                degree = degree,
                n_block = 3
            )
        }
        if (demean) {
            ## remove regional mean (uniform activity level) from each trial (all data)?
            B_prep <- lapply(B_prep, scale, center = TRUE, scale = FALSE)
        }
        if (divnorm) {
            ## divisive normalize each vertex by within-condition*run*session stdev (all data)?
            B_prep <- lapply(B_prep, function(betas){
                betas_i <- betas[, colnames(betas) %in% stimuli_eachrun]
                eps <- resid(.lm.fit(x = indicator(colnames(betas_i)), y = t(betas_i)))
                crossprod(diag(1/sqrt(Var(eps, 2))), betas)
            })
        }

        ## extract training dataset and implement initial subsampling:
        B_train <- B_prep[c("baseline_run1", "baseline_run2")]
        names(B_train) <- c("run1", "run2")
        b <- enlist(names(stimsets))
        for (stimset in names(stimsets)) {     
            b[[stimset]] <- enlist(variables)
            for (variable in variables) {
                ## skip bad combo:
                if (variable %in% c("target", "distractor") && stimset == "biaspc50") next                
                bi <- enlist(runs)
                for (run in runs) {
                    if (variable %in% c("target", "distractor")) to_drop <- stimuli_drop[[run]]
                    else if (variable == "congruency") to_drop <- ""
                    should_keep <- colnames(B_train[[run]]) %in% setdiff(stimsets[[stimset]], to_drop)
                    bii <- B_train[[run]][, should_keep]
                    colnames(bii) <- paste0(colnames(bii), "_", run)
                    bi[[run]] <- bii
                }
                b[[stimset]][[variable]] <- abind(bi, along = 2)  ## concatenate across runs
            }
        }

        ## resample trial indices for training models:
        if (roi_i == 1) {  ## same order for all ROIs: 
            idx <- enlist(names(stimsets))
            for (stimset in names(stimsets)) {  
                idx[[stimset]] <- enlist(variables)
                for (variable in variables) {

                    if (stimset == "biaspc50" & variable %in% c("target", "distractor")) next  ## skip bad combo
                    trial_labels <- colnames(b[[stimset]][[variable]])
                    set.seed(0)  ## ensure that target and distractor training sets are identical (b/c they can be)
                    if (variable %in% c("target", "distractor")) {
                        idx[[stimset]][[variable]] <- resample_idx(trial_labels, n_resamples, resample_to = 3)
                    } else if (variable == "congruency") {
                        if (stimset == "pc50") {
                            mat <- rbind(seq_along(trial_labels))  ## no resampling necessary
                            colnames(mat) <- trial_labels
                            idx[[stimset]][[variable]] <- mat
                        } else if (stimset %in% c("biaspc50", "bias")) {
                            ## conditions to downsample to 6
                            trial_idx6 <- which(trial_labels %in% stimuli_to6)
                            trial_labels_to6 <- trial_labels[trial_idx6]
                            to6 <- resample_idx(trial_labels_to6, n_resamples, resample_to = 6)  ## inds wrt trial_labels_to6
                            idx6 <- matrix(NA, ncol  = ncol(to6), nrow = nrow(to6))  ## but need to be wrt trial_inds...
                            dimnames(idx6) <- dimnames(to6)  ## so build new matrix idx6 and fill out
                            for (i in seq_along(trial_idx6)) idx6[to6 == i] <- trial_idx6[i]
                            ## conditions to downsample to 3
                            trial_idx3 <- which(trial_labels %in% stimuli_to3)
                            trial_labels_to3 <- trial_labels[trial_idx3]
                            to3 <- resample_idx(trial_labels_to3, n_resamples, resample_to = 3)  ## inds wrt trial_labels_to6
                            idx3 <- matrix(NA, ncol  = ncol(to3), nrow = nrow(to3))  ## but need to be wrt trial_inds...
                            dimnames(idx3) <- dimnames(to3)  ## so build new matrix idx6 and fill out
                            for (i in seq_along(trial_idx3)) idx3[to3 == i] <- trial_idx3[i]
                            ## rest of the items, keep all observations (no resamling necessary)
                            trial_idxrest <- which(!trial_labels %in% c(stimuli_to3, stimuli_to6))  ## all pc50, bias incon
                            idx_rest <- matrix(rep(trial_idxrest, n_resamples), nrow = n_resamples, byrow = TRUE)
                            colnames(idx_rest) <- trial_labels[trial_idxrest]
                            ## finally, store:
                            idx[[stimset]][[variable]] <- cbind(idx6, idx3, idx_rest)
                        }    
                    }

                }
            }
        }

        ## get regional means
        means <- lapply(B, colMeans)

        ## fit models and get predictions

        preds <- enlist(combo_paste(test_sessions, "_", variables, "_", names(stimsets)))
        for (test_session in test_sessions) {
            for (variable in variables) {
                for (stimset in names(stimsets)) {
                    
                    if (stimset == "biaspc50" & variable %in% c("target", "distractor")) next  ## skip bad combo
                    
                    nm <- paste0(test_session, "_", variable, "_", stimset)
                    
                    ## get train data:
                    idx_i <- idx[[stimset]][[variable]]
                    traindata <- b[[stimset]][[variable]]
                    colnames(traindata) <- gsub("_run[1-2]", "", colnames(traindata))
                    colnames(idx_i) <- gsub("_run[1-2]", "", colnames(idx_i))
                    y_train <- as.factor(get_variable(colnames(idx_i), variable))
                    
                    ## get test data:
                    testdata <- cbind(
                        B_prep[[paste0(test_session, "_run1")]], 
                        B_prep[[paste0(test_session, "_run2")]]
                    )
                    if (variable == "congruency") {
                        test_idx <- which(colnames(testdata) %in% stimsets$biaspc50)
                    }  else{
                        test_idx <- which(colnames(testdata) %in% stimsets[[stimset]])
                    }
                    testdata <- testdata[, test_idx]
                    y_test <- as.factor(get_variable(colnames(testdata), variable))

                    ## fit and predict:
                    n_resamples_i <- nrow(idx_i)  ## some don't require resampling
                    preds_i <- vector("list", n_resamples_i)
                    for (i in seq_len(n_resamples_i)) {

                        traindata_i <- traindata[, idx_i[i, ]]
                        stopifnot(colnames(traindata_i) == colnames(idx_i))

                        if (decoder == "pda") {
                            stopifnot("not yet configured for pda")
                        } else if (decoder == "svm") {
                            fit <- svm(x = t(traindata_i), y = y_train, kernel = "linear", scale = FALSE)
                        }

                        preds_i[[i]] <- summarize_predictions(fit = fit, newdata = t(testdata), decoder = decoder)
                        
                    }
                    preds_sum <- as.matrix(Reduce("+", preds_i) / n_resamples_i)
                    stopifnot(length(test_idx) == nrow(preds_sum))
                    
                    ## extract regional means:
                    mu <- c(means[[paste0(test_session, "_run1")]], means[[paste0(test_session, "_run2")]])
                    
                    preds[[nm]] <- data.table(
                        wave = unique(inputs$wave),
                        session = test_session,
                        variable = variable,
                        train_stimset = stimset,
                        roi = rois[roi_i],
                        trial = test_idx,
                        proj = preds_sum,
                        mu = mu[test_idx],
                        stimulus = rownames(preds_sum)
                    )
                    
                    stopifnot(length(test_idx) == nrow(preds[[nm]]))

                }
            }
        }  ## end fitting/prediction loop (session)

        preds <- preds[lengths(preds) > 0]
        res[[roi_i]] <- rbindlist(preds, fill = TRUE)
        B

        print(roi_i)

    }  ## end ROI loop

    res_dt <- rbindlist(res, id = "roi")
    # fname <- paste0(decoder, "__", atlas_name, "__", roi_col, "__nresamp", n_resamples, "__", glmname, "_", id, ".txt")
    # fwrite(res_dt, here("out", "res_decoding", fname))
    res_dt

}
stopCluster(cl)
time_end <- Sys.time()
print(time_end - time_beg)

outd <- rbindlist(out)
outd <- separate(outd, "roi", c("subj", "roi"), sep = "__")


## flags:
suffix_divnorm <- switch(divnorm + 1, "", "__divnormed")
suffix_demean <- switch(demean + 1, "", "__demeaned")
suffix_center <- switch(center + 1, "", "__centered")
suffix_detrend <- switch(detrend + 1, "", paste0("__detrended", degree))
fname <- here(
    "out", "decoding",
    paste0(
        decoder, "__", atlas_name, "__", roi_col, "__nresamp", n_resamples, "__", glmname, 
        suffix_divnorm, suffix_demean, suffix_center, suffix_detrend, 
        ".txt"
    )
)
fwrite(outd, fname)


## checking resampling code ----

# table(s(idx$bias))
# table(colnames(idx$pc50))
# table(colnames(idx$biaspc50))

# table(gsub("[A-z]", "", colnames(idx$bias)))
# table(gsub("[A-z]", "", colnames(idx$pc50)))
# table(gsub("[A-z]", "", colnames(idx$biaspc50)))

# table(get_variable(gsub("_run[0-9]", "", colnames(idx$bias)), "target"))
# table(get_variable(gsub("_run[0-9]", "", colnames(idx$pc50)), "target"))
# table(get_variable(gsub("_run[0-9]", "", colnames(idx$biaspc50)), "target"))

# table(get_variable(gsub("_run[0-9]", "", colnames(idx$bias)), "distractor"))
# table(get_variable(gsub("_run[0-9]", "", colnames(idx$pc50)), "distractor"))
# table(get_variable(gsub("_run[0-9]", "", colnames(idx$biaspc50)), "distractor"))

# table(get_variable(gsub("_run[0-9]", "", colnames(idx$bias)), "congruency"))
# table(get_variable(gsub("_run[0-9]", "", colnames(idx$pc50)), "congruency"))
# table(get_variable(gsub("_run[0-9]", "", colnames(idx$biaspc50)), "congruency"))
