#!/usr/bin/env bash Rscript --vanilla

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
#library(mda)
library(e1071)

## set variables

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

if (interactive()) {  ## add variables (potentially unique to this script) useful for dev
    glmname <- "lsall_1rpm"
    atlas_name <- "schaefer2018_17_200"
    #atlas_name <- "glasser2016"
    space <- "fsaverage5"
    roi_col <- "parcel"
    subjlist <- "wave1_unrel"
    subjects <- fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]]
    waves <- "wave1"
    n_cores <- 24
    n_resamples <- 1E2
    decoder <- "svm"
    # test_session <- "proactive"
} else {
    source(here("src", "parse_args.r"))
    print(args)
}

sessions <- c(train_session, test_sessions)
atlas <- load_atlas(atlas_name, space)
rois <- unique(atlas$key[[roi_col]])
## remove hippocampi from glasser atlas (not represented in fsaverages???)
if (atlas_name == "glasser2016" && grepl("fsaverage", space)) {
    rois <- rois[rois != "L_H" & rois != "R_H"]
}
roiset <- paste0(atlas_name, "_", roi_col)


## script-specific functions


summarize_predictions <- function(fit, newdata, labs, type = "contrast") {
    
    stopifnot(length(labs) == nrow(newdata))

    if (type == "contrast") {
        
        x <- attr(predict(fit, newdata = newdata, decision.values = TRUE), "decision.values")
        contrs <- colnames(x)
        l <- strsplit(colnames(x), "/")
        pos <- unlist(lapply(l, "[", 1))
        neg <- unlist(lapply(l, "[", 2))
        levs <- unique(labs)
        y <- matrix(NA, ncol = ncol(x), nrow = nrow(x), dimnames = dimnames(x))
        for (i in seq_along(contrs)) {
            pos_idx <- grep(pos[i], labs)
            neg_idx <- grep(neg[i], labs)
            y[pos_idx, i] <- x[pos_idx, i] - mean(x[neg_idx, i])
            y[neg_idx, i] <- mean(x[pos_idx, i]) - x[neg_idx, i]
        }
        ysum <- rowMeans(y, na.rm = TRUE)

    } else if (type == "probability") {
        
        x <- attr(predict(fit, newdata = newdata, probability = TRUE), "probabilities")
        x <- qlogis(x)
        ysum <- numeric(nrow(x))
        for (i in seq_len(nrow(x))) {
            idx <- grep(labs[i], colnames(x))
            ysum[i] <- mean(x[i, -idx]) - x[i, idx]
        }

    }
    
    ysum

}



## run ----


input <- construct_filenames_h5(
    prefix = "coefs", subjects = subjects, waves = waves, sessions = sessions, rois = rois, runs = runs, 
    glmname = glmname, prewh = "none"
)
## identify and remove subjects with missing data:
input[, file_exists := file.exists(file_name)]
subjs_missing <- input[file_exists == FALSE, unique(subj)]
print(noquote(paste0("removing subjects due to missing data: ", subjs_missing)))
input <- input[!subj %in% subjs_missing]
## make larger chunks for each core by redefining looping variable (have each core loop over ROIs):
## each row of input_subset defines the location of data for a single subject*wave*session*roi (both runs):
input[, g := subj]  ## run subjects in parallel (ROI in serial)
g <- unique(input$g)


time_beg <- Sys.time()
cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
out <- foreach(ii = seq_along(g), .verbose = TRUE, .inorder = FALSE) %dopar% {

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
        
        ## extract training dataset and implement initial subsampling:

        B_train <- B[c("baseline_run1", "baseline_run2")]
        names(B_train) <- c("run1", "run2")
        b <- enlist(names(stimsets))
        for (stimset in names(stimsets)) { 
            ## extract baseline session data and drop conditions that cannot be resampled :(
            bi <- enlist(runs)
            for (run in runs) {
                to_drop <- switch(stimset, bias = stimuli_drop[[run]], pc50 = stimuli_drop[[run]], biaspc50 = "")
                should_keep <- colnames(B_train[[run]]) %in% setdiff(stimsets[[stimset]], to_drop)
                bii <- B_train[[run]][, should_keep]
                colnames(bii) <- paste0(colnames(bii), "_", run)
                bi[[run]] <- bii
            }
            b[[stimset]] <- abind(bi, along = 2)  ## concatenate across runs
        }

        ## resample trial indices for training models:

        if (roi_i == 1) {  ## same order for all ROIs: 
            idx <- enlist(names(stimsets))
            for (stimset in names(stimsets)) {  
                trial_labels <- colnames(b[[stimset]])
                set.seed(0)
                
                if (stimset %in% c("bias", "pc50")) {  ## target and distractor models

                    idx[[stimset]] <- resample_idx(trial_labels, n_resamples, resample_to = 3)

                } else if (stimset == "biaspc50") {  ## congruency model
                    
                    trial_idx6 <- which(trial_labels %in% stimuli_to6)
                    trial_labels_to6 <- trial_labels[trial_idx6]
                    to6 <- resample_idx(trial_labels_to6, n_resamples, resample_to = 6)  ## inds wrt trial_labels_to6
                    idx6 <- matrix(NA, ncol  = ncol(to6), nrow = nrow(to6))  ## but need to be wrt trial_inds...
                    dimnames(idx6) <- dimnames(to6)  ## so build new matrix idx6 and fill out
                    for (i in seq_along(trial_idx6)) idx6[to6 == i] <- trial_idx6[i]
                    
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
                    
                    idx[[stimset]] <- cbind(idx6, idx3, idx_rest)

                }
            }
        }

        ## scale by mean of mean lengths?
        ## scaling condition-patterns separately can drive dependence on mean pattern...
        d <- enlist(test_sessions)
        for (test_session in test_sessions) {
            nms <- paste0(test_session, "_", runs)
            b_test <- B[nms]
            b_test <- abind(b_test, along = 2)  ## concatenate across runs
            trial_sd <- sqrt(Var(b_test[, colnames(b_test) %in% stimsets$pc50], 2))
            trial_congruency <- get_variable(names(trial_sd), "congruency")
            sd_bar <- mean(vapply(split(trial_sd, trial_congruency), mean, numeric(1)))            
            d[[test_session]] <- b_test / sd_bar
        }


        ## fit models and get predictions

        preds <- enlist(combo_paste(test_sessions, "_", variables, "_", names(stimsets)))
        for (test_session in test_sessions) {
            for (variable in variables) {
                for (stimset in names(stimsets)) {
                    
                    nofit <- 
                        (variable %in% c("target", "distractor") & stimset == "biaspc50") |
                        (variable == "congruency" & stimset %in% c("bias", "pc50")) |
                        (test_session == "reactive" & stimset == "bias")
                    if (nofit) next

                    nm <- paste0(test_session, "_", variable, "_", stimset)
                    ## get train data:
                    idx_i <- idx[[stimset]]
                    traindata <- b[[stimset]]
                    colnames(traindata) <- gsub("_run[1-2]", "", colnames(traindata))
                    colnames(idx_i) <- gsub("_run[1-2]", "", colnames(idx_i))
                    y_train <- as.factor(get_variable(colnames(idx_i), variable))
                    ## get test data:
                    testdata <- d[[test_session]]
                    use_these <- switch(variable, 
                        target = stimsets[[stimset]], distractor = stimsets[[stimset]], congruency = stimsets$pc50)
                    testdata <- testdata[, colnames(testdata) %in% use_these]
                    y_test <- as.factor(get_variable(colnames(testdata), variable))

                    preds_i <- vector("list", n_resamples)
                    for (i in seq_len(n_resamples)) {

                        traindata_i <- traindata[, idx_i[i, ]]
                        stopifnot(colnames(traindata_i) == colnames(idx_i))                        
                        
                        ## optional: aggregate into pseudotrials? (for no aggregation, set n_tpp = 1)
                        #d_train_i <- sample_pseudotrials(d_train_i, y_train, n_tpp)
                        
                        ## fit model
                        
                        if (decoder == "pda") {

                        } else if (decoder == "svm") {
                            fit <- svm(x = t(traindata_i), y = y_train, kernel = "linear", probability = TRUE, scale = FALSE)
                        }

                        ## get predictions
                            
                        if (decoder == "pda") {

                        } else if (decoder == "svm") {
                            ## function: fit, newdata, y_test
                            distn <- summarize_predictions(fit, t(testdata), y_test, "contrast")
                            #logit <- summarize_predictions(fit, t(testdata), y_test, "probability")
                            #preds_i[[i]] <- cbind(distn, logit)
                            preds_i[[i]] <- distn
                        }
                        
                    }  ## end resamples

                    preds[[nm]] <- Reduce("+", preds_i) / n_resamples

                }
            }
        }  ## end fitting/prediction loop (session)

        preds <- preds[lengths(preds) > 0]
        res[[roi_i]] <- rbindlist(lapply(preds, as.data.table), id = "session_variable_stimset")

        #print(roi_i)

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
outd <- separate(outd, "session_variable_stimset", c("session", "variable", "stimset"), sep = "_")
outd <- rename(outd, distn = V1)
fname <- paste0(decoder, "__", atlas_name, "__", roi_col, "__nresamp", n_resamples, "__", glmname, ".txt")
fwrite(outd, here("out", "res_decoding", fname))


## checking resampling code ----

# table(colnames(idx$bias))
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
