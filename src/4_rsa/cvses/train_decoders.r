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
library(mda)
library(e1071)


## set variables

task <- "Stroop"
variables <- c("target", "distractor", "congruency")
stimsets <- c("bias", "pc50", "biaspc50")
train_session <- "baseline"
test_sessions <- c("proactive")
# stimuli <- list(
#     pc50 = combo_paste(colors_pc50, words_pc50),
#     bias = combo_paste(colors_bias, words_bias)
# )
trialcounts <- fread(here("in", "trialcounts.csv"))
## stimuli to drop from training set:
stimuli_drop_balance <- list(
    run1 = c("blackBLACK", "pinkPINK", "redRED", "whiteWHITE"),  ## items to drop in run 1 (and downsample in run 2)
    run2 = c("yellowYELLOW", "greenGREEN", "blueBLUE", "purplePURPLE")  ## items to drop in run 2 (downsample run 1)
)
## stimuli to drop when centering target and distractor models of PC50 data:
stimuli_drop_centering <- list(
    run1 = trialcounts[pc == "pc50" & count == 0 & run == "run1"]$stimulus,
    run2 = trialcounts[pc == "pc50" & count == 0 & run == "run2"]$stimulus
)
stimsets <- list(
    bias = trialcounts[pc == "bias", unique(stimulus)],
    pc50 = trialcounts[pc == "pc50", unique(stimulus)]
)


#dsets <- c("bias", "pc50", "biaspc50")
# models <- data.frame(
#     variable = c("target", "distractor", "target", "distractor", "congruency"),
#     stimset = c("bias", "pc50", "bias", "pc50", "bias+pc50")
# )

# trialcounts[pc == "pc50", sum(count), by = c("stimulus", "session")] %>% dcast(stimulus ~ session)
# trialcounts[pc == "pc50", sum(count), by = c("run", "session")] %>% dcast(run ~ session)
# stimuli_drop
# c(, c("blackBLACK", "pinkPINK", "redRED", "whiteWHITE"))


# trialcounts <- melt(trialcounts, id.vars = c("stimulus", "pc"), value.name = "count", variable.name = "session")
# trialcounts$run <- paste0("run", gsub("[a-z]", "", trialcounts$session))
# trialcounts$session <- gsub("[0-9]", "", trialcounts$session)
# fwrite(trialcounts, here("in", "trialcounts.csv"))


if (interactive()) {  ## add variables (potentially unique to this script) useful for dev
    glmname <- "lsall_1rpm"
    atlas_name <- "glasser2016"
    space <- "fsaverage5"
    roi_col <- "parcel"
    subjlist <- "wave1"
    subjects <- fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]]
    waves <- "wave1"
    roi_i <- 1
    ii <- 1
    n_cores <- 20
    n_resamples <- 1E2
    decoder <- "pda"
    shrinkage_factor <- 1
    n_tpp <- 2  ## n trials per pseudotrial
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

get_variable <- function(x, variable) {
    if (variable == "target") {
        res <- gsub("[[:upper:]]", "", x)
    } else if (variable == "distractor") {
        res <- gsub("[[:lower:]]", "", x)
    } else if (variable == "congruency") {
        is_congr <- gsub("[[:upper:]]", "", x) == tolower(gsub("[[:lower:]]", "", x))
        res <- ifelse(is_congr, "congr", "incon")
    } else res <- "variable must be target, distractor, or congruency"

    res

}

# preproc <- function(l, center = TRUE, scale = TRUE, return_matrix = TRUE) {
#     stopifnot(is.list(l) && length(l) == 2)
    
#     if (center) {

#         ## average over t*d in each run:
#         mu_td_run <- lapply(l, function(x) average(x, colnames(x)))
#         ## concatenate t*d*run vectors columnwise into single matrix:
#         mu_td_run <- abind(mu_td_run, along = 2)
#         ## average over run by t*d (for t*d in both runs, i.e., congruents):
#         mu_td <- average(mu_td_run, colnames(mu_td_run))
#         mu_bar <- rowMeans(mu_td)  ## now over t*d
#         l <- lapply(l, "-", mu_bar)

#     }

#     if (scale) {
#         dn <- lapply(l, dimnames)
#         l <- lapply(l, scale2unit)
#         for (i in seq_along(l)) dimnames(l[[i]]) <- dn[[i]]
#     }

#     if (return_matrix) l <- abind(l, along = 2)
    
#     l

# }

# resample_idx_congr <- function(trial_labels_bias, n_resamples) {

#     ## only bias congruent items need to be resampled.
#     ## blueBLUE and purplePURPLE occur 15 times in run 1 and 12 times in run 2. 
#     ## they need to be resampled to 6 times in run 1 and 3 times in run 2.
#     ## redRED and whiteWHITE occur 15 times in run 2 and 12 times in run 1. 
#     ## they need to be resampled to 6 times in run 2 and 3 times in run 1.
#     ## all other items (pc50, bias incongruent) do not need to be resampled.
#     ## this stratified resampling ensures that each target or distractor occurs equally often as congruent
#     ## and incongruent, and that congruent vs incongruent trials occur with equal frequency.

#     resamp_info <- data.frame(
#         item = rep(c("blueBLUE", "purplePURPLE", "redRED", "whiteWHITE"), 2),
#         run = c(rep("run1", 4), rep("run2", 4)),
#         to = c(6, 6, 3, 3, 3, 3, 6, 6),
#         stringsAsFactors = FALSE
#     )
#     idx_resamp <- enlist(paste0(resamp_info$item, "_", resamp_info$run))
#     for (i in seq_len(nrow(resamp_info))) {
#         item <- resamp_info$item[i]
#         run <-  resamp_info$run[i]
#         to <- resamp_info$to[i]
#         labs <- trial_labels_bias[[paste0(train_session, "_", run)]]
#         idxmat <- t(replicate(n_resamples, resample(which(item == labs), to)))
#         colnames(idxmat) <- rep(item, ncol(idxmat))
#         idx_resamp[[i]] <- idxmat
#     }
#     ## now concatenate different resampled items (within run):
#     idx_bias <- list(
#         run1 = abind(idx_resamp[resamp_info$run == "run1"], along = 2),
#         run2 = abind(idx_resamp[resamp_info$run == "run2"], along = 2)
#     )
#     ## and concatenate all other items:
#     for (run_i in 1:2) {
#         is_allelse <- !trial_labels_bias[[run_i]] %in% resamp_info$item
#         idx_allelse <- matrix(rep(which(is_allelse), n_resamples), ncol = sum(is_allelse), byrow = TRUE)
#         colnames(idx_allelse) <- trial_labels_bias[[run_i]][is_allelse]
#         idx_bias[[run_i]] <- cbind(idx_bias[[run_i]], idx_allelse)
#     }

#     idx_bias

# }

summarize_predictions <- function(x, labs, into = "confusion") {
    labs <- as.character(labs)
    if (into == "confusion") {
        A <- averaging_matrix(labs)
        res <- crossprod(A, x)
        names(dimnames(res)) <- c("true", "predicted")
    } else if (into == "contrast") {
        res <- numeric(nrow(x))
        for (i in seq_len(nrow(x))) {
            idx <- grep(labs[i], colnames(x))
            res[i] <- mean(x[i, -idx]) - x[i, idx]
        }
    }
    res
}

# pdist2 <- function(A,B) {
  
#   ## this function works on matrices A and B.
#   ## A and B should be matrices with observations as rows and features as columns.
#   ## this function computes the squared euclidean distances between each row of A and each row of B.
#   ## the output is therefore a matrix of size nrow(A)*nrow(B).
  
#   ## this is an efficient implementation of the pdist::pdist function.
#   ## see links for more information:
#   ## https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
#   ## https://blog.smola.org/post/969195661/in-praise-of-the-second-binomial-formula
  
  
#   an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
#   bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  
#   m = nrow(A)
#   n = nrow(B)
  
#   tmp = matrix(rep(an, n), nrow=m) 
#   tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  
#   tmp - 2 * tcrossprod(A,B)  ## squared euclidean distance
  
# }

# mahals <- function(fit, newdata) {
#     variates <- predict(fit, newdata, type = "variates")
#     pdist2(variates, fit$means)
# }


# sample_pseudotrials <- function(dat, labs, n_tpp) {
    
#     idx <- split(seq(ncol(dat)), labs)
#     idx_shuffle <- enlist(names(idx))
#     labs_pseudo <- enlist(names(idx))
#     for (cond_i in seq_along(idx)) {
#         idx_cond <- idx[[cond_i]]
#         idx_shuffle[[cond_i]] <- sample(idx_cond)
#         labs_pseudo[[cond_i]] <- paste0(names(idx)[cond_i], ceiling(seq_along(idx_cond)/n_tpp))
#     }
#     idx_shuffle <- unlist(idx_shuffle, use.names = FALSE)
#     labs_pseudo <- unlist(labs_pseudo, use.names = FALSE)
#     dat <- average(dat[, idx_shuffle], labs_pseudo)
#     colnames(dat) <- gsub("[0-9]", "", colnames(dat))

#     dat

# }

















stimuli_train <- list(
    bias = list(
        run1 = c("blueBLUE", "blueWHITE", "purpleBLUE", "purplePURPLE", "redBLUE", "redPURPLE", "whitePURPLE", "whiteRED"),
        run2 = c("bluePURPLE", "blueRED", "purpleRED", "purpleWHITE", "redRED", "redWHITE", "whiteBLUE", "whiteWHITE")
    ),
    pc50 = list(
        run1 = c("blackPINK", "blackYELLOW", "yellowGREEN", "yellowYELLOW", "greenBLACK", "greenGREEN", "pinkGREEN", "pinkYELLOW"),
        run2 = c("blackBLACK","blackGREEN","yellowBLACK","yellowPINK","greenPINK","greenYELLOW","pinkBLACK","pinkPINK")
    ),
    biaspc50 = list(
        run1 = c(
            "blackPINK", "blackYELLOW", "yellowGREEN", "yellowYELLOW", "greenBLACK", "greenGREEN", "pinkGREEN", "pinkYELLOW",
            "blueBLUE", "blueWHITE", "purpleBLUE", "purplePURPLE", "redBLUE", "redPURPLE", "whitePURPLE", "whiteRED"
            ),
        run2 = c(
            "blackBLACK","blackGREEN","yellowBLACK","yellowPINK","greenPINK","greenYELLOW","pinkBLACK","pinkPINK",
            "bluePURPLE", "blueRED", "purpleRED", "purpleWHITE", "redRED", "redWHITE", "whiteBLUE", "whiteWHITE"
            )
    )
)

center_training_data <- function(b, session, stimset) {
    
}


## .center4targetdistractor
## .center4conruency_model

center_data <- function(b, session, stimset, .stimuli_train = stimuli_train[[stimset]], .runs = runs) {

    ## given a session, stimset, and model
    ## if stimset == "biaspc50", then variable %in% c("target", "distractor")
    #session <- "baseline"
    #stimset <- "pc50"
    #variable <- "target"
    ## make explicit list of stimuli*runs to use per dset (pc50, bias, pc50+bias)
    
    stopifnot(is.list(b) && length(b) == 2 && identical(names(b), paste0(session, "_", .runs)))

    ## extract trialtypes to use:
    for (run in .runs) {
        nm <- paste0(session, "_", run)
        bi <- b[[nm]]
        stim_keep <- .stimuli_train[[run]]
        b[[nm]] <- bi[, colnames(bi) %in% stimuli_keep_i]
    }
    b <- abind(b, along = 2)  ## concatenate across runs
    
    ## get dataset
    if (session == "baseline") {
        dset <- b  ## baseline ("training") uses only subsampled data
    } else {
        dset <- abind(B[b_idx], along = 0)  ## proactive and reactive ('test') use 
    }

    ## compute mean of means (via downsampling with baseline/proactive, via all data with reactive)
    ## and center dataset
    if (session %in% c("baseline", "proactive")) {
        to3 <- resample_idx(colnames(b), n_resamples, resample_to = 3)
        to1 <- resample_idx(colnames(b), n_resamples, resample_to = 1)
        to3 <- to3[, get_variable(colnames(to3), "congruency") == "congr"]
        to1 <- to1[, get_variable(colnames(to1), "congruency") == "incon"]
        mu_idx <- cbind(to3, to1)
        mu <- vector("list", n_resamples)
        for (i in seq_len(n_resamples)) {
            bi <- b[, mu_idx[i, ]]
            mu[[i]] <- rowMeans(average(bi, colnames(bi)))
        }
        dset_c <- lapply(mu, function(x) dset - x)  ## center each condition at mean of means

    } else if (session == "reactive") {
        mu <- rowMeans(average(bi, colnames(bi)))
        dset_c <- dset - mu
    }  ## if baseline/proactive, dset_c is list of centered matrices; if reactive, dset_c is matrix

    dset_c
    
}



## execute ----

input <- construct_filenames_h5(
    prefix = "coefs", subjects = subjects, waves = waves, sessions = sessions, rois = rois, runs = runs, 
    glmname = glmname, prewh = "none"
)
## make larger chunks for each core by redefining looping variable (have each core loop over ROIs):
## each row of input_subset defines the location of data for a single subject*wave*session*roi (both runs):
input[, g := subj]  ## run subjects in parallel (ROI in serial)
g <- unique(input$g)

time_beg <- Sys.time()
cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
res <- foreach(ii = seq_along(g) %dopar% {

    id <- g[ii]
    inputs <- input[g == id]
    #ses <- unique(inputs$session)

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
        
        ## preprocess data
        
        ## drop conditions that cannot be resampled :(
        #train_idx <- grep(train_session, names(B))
        #lapply(B[train_idx], function(x) x[, colnames(x) %in% stimuli[[stimset]]])  
        
        for (stimset in names(stimsets)) {

            b_idx <- grep(session, names(B))
            b <- B[b_idx]
            bc <- center_data(b, session, stimset)

        }
        
        ## resample trial indices for training models:
        if (roi_i == 1) {  ## same order for all ROIs: 
            
            trial_labels <- colnames(d[[train_session]])
            ## get trial indices:    
            #trial_labels_pc50 <- lapply(d[[train_session]]$pc50, colnames)
            #trial_labels_bias <- lapply(d[[train_session]]$bias, colnames)
            for (stimset in stimsets) {
                set.seed(0)

                resample_idx(trial_labels, )
                
                train_idx$td <- list(
                    pc50 = lapply(trial_labels, resample_idx, n_resamples = n_resamples, resample_to = 3),
                    bias = lapply(trial_labels_bias, resample_idx, n_resamples = n_resamples, resample_to = 3)
                )
            train_idx$c <- resample_idx_congr(trial_labels_bias, n_resamples)  ## no resampling necessary for c_pc50


            }













            #train_idx <- enlist(c("td", "c"))  ## targetdistractor or congruency
            ## td (variable %in% target, distractor):
            ## 
            set.seed(0)
            train_idx$td <- list(
                pc50 = lapply(trial_labels_pc50, resample_idx, n_resamples = n_resamples, resample_to = 3),
                bias = lapply(trial_labels_bias, resample_idx, n_resamples = n_resamples, resample_to = 3)
            )
            train_idx$c <- resample_idx_congr(trial_labels_bias, n_resamples)  ## no resampling necessary for c_pc50

        }  

        for (model_i in seq_len(nrow(models))) {

            variable <- models$variable[model_i]
            stimset <- models$stimset[model_i]

            if (stimset %in% c("pc50", "bias")) {
                d_train <- d[[train_session]][[stimset]]
                train_idx_i <- train_idx$td[[stimset]]
                d_test <- d$proactive[[stimset]]
                names(d_test) <- c("run1", "run2")
                d_test_i <- lapply(d_test, function(x) x[, colnames(x) %in% stimuli[[stimset]]])
            } else if (stimset == "bias+pc50") {
                #d_train <- 
            }

            #fits <- vector("list", n_resamples)
            for (i in seq_len(n_resamples)) {
                
                ## train model

                run1_idx <- train_idx_i[[1]][i, ]
                run2_idx <- train_idx_i[[2]][i, ]
                d_train_i <- cbind(d_train[[1]][run1_idx, ], d_train[[2]][run2_idx, ])
                y_train <- get_variable(colnames(d_train_i), variable)

                ## optional: aggregate into pseudotrials? (for no aggregation, set n_tpp = 1)
                d_train_i <- sample_pseudotrials(d_train_i, y_train, n_tpp)
                
                if (decoder == "pda") {
                    fit <- fda(as.factor(colnames(d_train_i)) ~ t(scale2unit(d_train_i)), method = gen.ridge, lambda = shrinkage_factor)
                } else if (decoder == "svm") {
                    fit <- svm(x = t(scale2unit(training_set_i_t)), y = as.factor(colnames(d_train_i)), kernel = "linear", probability = TRUE)
                }          
                
                d_test_i <- preproc(d_test)  ## center at mean of means
                y_test <- get_variable(colnames(d_test_i), variable)
                d_test_i <- sample_pseudotrials(d_test_i, y_test, n_tpp)
                
     

            }


        }



        }



        
        res_stimset <- enlist(names(stimuli))
        for (stimset_i in seq_along(stimuli)) {

            stimset <- names(stimuli)[stimset_i]
            
            ## train model (baseline data) ----
            
            B_train_i <- center_mom(B_train[[stimset]])    ## center at mean of means
            
            ## fit model:
            fits_t <- vector("list", n_resamples)
            fits_d <- vector("list", n_resamples)
            for (i in seq_len(n_resamples)) {

                training_set_i <- cbind(
                    B_train_i$run1[, downsamp_idx[[stimset]]$run1[i, ]], 
                    B_train_i$run2[, downsamp_idx[[stimset]]$run2[i, ]]
                    )

                ## optional: aggregate into pseudotrials? (for no aggregation, set n_tpp = 1)
                training_set_i_d <- sample_pseudotrials(training_set_i, get_upper(colnames(training_set_i)), n_tpp)
                training_set_i_t <- sample_pseudotrials(training_set_i, get_lower(colnames(training_set_i)), n_tpp)

                ## scale (pseudo-)trial vectors to unit lengths:: TODO: change scale2unit to preserve dimnames!!
                # training_set_i_d <- scale2unit(training_set_i_d)
                # training_set_i_t <- scale2unit(training_set_i_t)

                ## train models
                if (decoder == "pda") {
                    fits_t[[i]] <- fda(as.factor(colnames(training_set_i_t)) ~ t(scale2unit(training_set_i_t)), method = gen.ridge, lambda = shrinkage_factor)
                    fits_d[[i]] <- fda(as.factor(colnames(training_set_i_d)) ~ t(scale2unit(training_set_i_d)), method = gen.ridge, lambda = shrinkage_factor)
                } else if (decoder == "svm") {
                    fits_t[[i]] <- svm(x = t(scale2unit(training_set_i_t)), y = as.factor(colnames(training_set_i_t)), kernel = "linear", probability = TRUE)
                    fits_d[[i]] <- svm(x = t(scale2unit(training_set_i_d)), y = as.factor(colnames(training_set_i_d)), kernel = "linear", probability = TRUE)
                }
                
            }

            ## test model (proactive, reactive data) ----

            res_session <- enlist(test_sessions)
            for (test_session in test_sessions) {
        
                ## prepare testing data
                B_test <- B[grep(test_session, names(B))]
                names(B_test) <- c("run1", "run2")
                B_test_i <- lapply(B_test, function(x) x[, colnames(x) %in% stimuli[[stimset_i]]])
                B_test_i <- center_mom(B_test_i)  ## center at mean of means
                B_test_i <- abind(B_test_i, along = 2)
                ## optional: aggregate into pseudotrials? (for no aggregation, set n_tpp = 1)
                B_test_d <- sample_pseudotrials(B_test_i, get_upper(colnames(B_test_i)), n_tpp)
                B_test_t <- sample_pseudotrials(B_test_i, get_lower(colnames(B_test_i)), n_tpp)
                
                ## compute predictions

                if (decoder == "pda") {
                    preds_t <- lapply(fits_t, mahals, newdata = t(scale2unit(B_test_t)))
                    preds_d <- lapply(fits_d, mahals, newdata = t(scale2unit(B_test_d)))
                } else if (decoder == "svm") {
                    #preds_t <- lapply(fits_t, function(x) attr(predict(x, newdata = t(scale2unit(B_test_t)), probability = TRUE), "probabilities"))
                    #preds_d <- lapply(fits_d, function(x) attr(predict(x, newdata = t(scale2unit(B_test_d)), probability = TRUE), "probabilities"))
                    ## TODO: if going to use, add in qlogis transform
                }
                
                ## summarize predictions

                ## aggregate predictions over resampled models:
                preds_bar <- list(
                    target = Reduce("+", preds_t) / n_resamples,
                    distractor = Reduce("+", preds_d) / n_resamples
                )                
                ## get testing labels
                y_test <- list(
                    target = colnames(B_test_t),
                    distractor = colnames(B_test_d)
                )

                #confusions <- Map(summarize_predictions, preds_bar, y_test, into = "confusion")
                contrs <- abind(Map(summarize_predictions, preds_bar, y_test, into = "contrast"), rev.along = 0)
                predclass <- abind(lapply(preds_bar, function(x) colnames(x)[apply(x, 1, which.min)]), rev.along = 0)
                true <- abind(y_test, rev.along = 0)
                out <- data.table(contrs, predclass, true)
                names(out) <- c(
                    "contr_target", "contr_distractor", "pred_target", "pred_distractor", "true_target", 
                    "true_distractor"
                    )
                res_session[[test_session]] <- out
            }
            
            res_stimset[[stimset_i]] <- rbindlist(res_session, id = "session")

        }

        res[[roi_i]] <- rbindlist(res_stimset, id = "stimset")

    }

    res <- rbindlist(res, id = "subj__roi")
    separate(res, col = "subj__roi", c("subj", "roi"), "__")  ## return

}
stopCluster(cl)
time_end <- Sys.time()
time_end - time_beg


d <- rbindlist(res)

#saveRDS(d, here("out", "pda_lambda=1_pseudo=1_nresamp=1000_stimset-sep.RDS"))
#saveRDS(d, here("out", "pda_lambda=1_pseudo=2_nresamp=100_stimset-sep.RDS"))
saveRDS(d, here("out", "pda_lambda=100_pseudo=2_nresamp=100_stimset-sep_meas-mah.RDS"))

d <- readRDS(here("out", "pda_lambda=100_pseudo=2_nresamp=100_stimset-sep_meas-mah.RDS"))

d_sum <- d %>%
    group_by(subj, roi, stimset) %>%
    summarize(
        cr_target = mean(pred_target == true_target),
        cr_distractor = mean(pred_distractor == true_distractor),
        contr_distractor = mean(contr_distractor),
        contr_target = mean(contr_target)
        )

d_sum <- d_sum %>%
    as.data.table %>%
    melt(
        id.vars = c("subj", "roi", "stimset"), 
        measure.vars = patterns(cr = "^cr", contr = "^contr"),
        variable.name = "model"
        )
d_sum$model <- c("target", "distractor")[d_sum$model]

d_sum2 <- d_sum %>%
    group_by(roi, stimset, model) %>%
    summarize(
        contr_p = t.test(contr, alternative = "greater")$p.value, 
        cr_p = t.test(cr - 0.25, alternative = "greater")$p.value,
        contr_t = t.test(contr)$statistic, 
        cr_t = t.test(cr - 0.25)$statistic,
        cr = mean(cr), 
        contr = mean(contr)
    ) %>%
    group_by(stimset, model) %>%
    mutate(
        contr_p_fdr = p.adjust(contr_p, "fdr"), 
        cr_p_fdr = p.adjust(cr_p, "fdr")
        )

d_sum2 %>% filter(contr_p_fdr < 0.05) %>% View
d_sum2 %>% filter(cr_p_fdr < 0.05) %>% View

