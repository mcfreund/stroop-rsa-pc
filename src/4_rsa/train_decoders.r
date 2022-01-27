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

get_upper <- function(x) gsub("[[:lower:]]", "", x)
get_lower <- function(x) gsub("[[:upper:]]", "", x)

center_mom <- function(l) {
    ## average over t*d in each run:
    mu_td_run <- lapply(l, function(x) average(x, colnames(x)))
    ## concatenate t*d*run vectors columnwise into single matrix:
    mu_td_run <- abind(mu_td_run, along = 2)
    ## average over run by t*d (for t*d in both runs, i.e., congruents):
    mu_td <- average(mu_td_run, colnames(mu_td_run))
    mu_bar <- rowMeans(mu_td)  ## now over t*d
    lapply(l, "-", mu_bar)
}

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

pdist2 <- function(A,B) {
  
  ## this function works on matrices A and B.
  ## A and B should be matrices with observations as rows and features as columns.
  ## this function computes the squared euclidean distances between each row of A and each row of B.
  ## the output is therefore a matrix of size nrow(A)*nrow(B).
  
  ## this is an efficient implementation of the pdist::pdist function.
  ## see links for more information:
  ## https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
  ## https://blog.smola.org/post/969195661/in-praise-of-the-second-binomial-formula
  
  
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  
  m = nrow(A)
  n = nrow(B)
  
  tmp = matrix(rep(an, n), nrow=m) 
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  
  tmp - 2 * tcrossprod(A,B)  ## squared euclidean distance
  
}

mahals <- function(fit, newdata) {
    variates <- predict(fit, newdata, type = "variates")
    pdist2(variates, fit$means)
}


sample_pseudotrials <- function(dat, labs, n_tpp) {
    
    idx <- split(seq(ncol(dat)), labs)
    idx_shuffle <- enlist(names(idx))
    labs_pseudo <- enlist(names(idx))
    for (cond_i in seq_along(idx)) {
        idx_cond <- idx[[cond_i]]
        idx_shuffle[[cond_i]] <- sample(idx_cond)
        labs_pseudo[[cond_i]] <- paste0(names(idx)[cond_i], ceiling(seq_along(idx_cond)/n_tpp))
    }
    idx_shuffle <- unlist(idx_shuffle, use.names = FALSE)
    labs_pseudo <- unlist(labs_pseudo, use.names = FALSE)
    dat <- average(dat[, idx_shuffle], labs_pseudo)
    colnames(dat) <- gsub("[0-9]", "", colnames(dat))

    dat

}

## set variables

G <- c(colors_bias, colors_pc50)
items_drop <- list(
    run1 = c("blackBLACK", "pinkPINK", "redRED", "whiteWHITE"),  ## items to drop in run 1 and downsample in run 2
    run2 = c("yellowYELLOW", "greenGREEN", "blueBLUE", "purplePURPLE")  ## items to drop in run 2, downsample run 1
)
task <- "Stroop"

if (interactive()) {  ## add variables (potentially unique to this script) useful for dev
    glmname <- "lsall_1rpm"
    atlas_name <- "glasser2016"
    space <- "fsaverage5"
    roi_col <- "parcel"
    subjlist <- "wave1"
    subjects <- fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]]
    waves <- "wave1"
    #sessions <- c("baseline", "proactive")
    measure <- "crcor"  ## crcor
    prewh <- "none"
    ttype_subset <- "all"
    ii <- 1
    roi_i <- 1
    i <- 1
    n_cores <- 20
    run_i <- 1
    n_resamples <- 1E2
    overwrite <- TRUE
    stimset_i = 1
    decoder <- "pda"
    test_session <- "proactive"
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

stimuli <- list(
    pc50 = combo_paste(colors_pc50, words_pc50),
    bias = combo_paste(colors_bias, words_bias)
)

train_session <- "baseline"
test_sessions <- c("proactive")
sessions <- c(train_session, test_sessions)

shrinkage_factor <- 1
n_tpp <- 2  ## n trials per pseudotrial

## execute ----


input <- construct_filenames_h5(
    prefix = "coefs", subjects = subjects, waves = waves, sessions = sessions, rois = rois, runs = runs, 
    glmname = glmname, prewh = prewh
)
input[, g := paste0(subj, "__", wave, "__", session, "__", roi)]
l <- split(input, by = "g")


## skip computation of existing output files if overwrite == FALSE:

# fname <- construct_filename_rdm(
#     measure = measure, glmname = glmname, ttype_subset = ttype_subset, roiset = roiset, prewh = prewh
#     )  ## output file

# if (!overwrite) {
#     if (file.exists(fname)) {
#         group_name <- combo_paste(subjects, waves, sessions, sep = "__")
#         group_exists <- group_name %in% h5ls(fname)$name
#         if (any(group_exists)) {
#             existing_group_idx <- grep(paste0(group_name[group_exists], collapse = "|"), names(l))
#             l <- l[-existing_group_idx]
#             if (length(l) == 0L) stop("All to-be-written files already exist. Exiting...")
#         }
#     }
# }



## make larger chunks for each core by redefining looping variable (have each core loop over ROIs):

input_subset <- rbindlist(l)
## each row of input_subset defines the location of data for a single subject*wave*session*roi (both runs):
input_subset[, g := subj]  ## run subjects in parallel (ROI in serial)
l <- split(input_subset, by = "g")

l <- l[1:20]


time_beg <- Sys.time()
cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
res <- foreach(ii = seq_along(l)) %dopar% {

    inputs <- l[[ii]]
    ses <- unique(inputs$session)
    id <- names(l)[ii]

    res <- enlist(paste0(id, "__", rois))
    for (roi_i in seq_along(rois)) {
        
        input_val <- inputs[roi == rois[roi_i]]

        ## read betas and subset vertices to common "good" set
        
        B <- enlist(combo_paste(sessions, "_", runs))
        is_bad_vert <- enlist(combo_paste(sessions, "_", runs))
        for (session_nm in sessions) {
            for (run_nm in runs) {
                sesrun <- paste0(session_nm, "_", run_nm)
                nms <- input_val[run == run_nm & session == session_nm, c("file_name", "dset_name")]
                B[[sesrun]] <- read_dset(nms$file_name, nms$dset_name)
                is_bad_vert[[sesrun]] <- is_equal(Var(B[[sesrun]], 1), 0)  ## verts with no BOLD
            }
        }
        is_good_vert <- !Reduce("|", is_bad_vert)
        if (sum(is_good_vert) == 0) {
            stop(c("no good verts: ", paste0(unique(input_val$subj), " ", unique(input_val$roi), sep = " ")))
        } else {
            B <- lapply(B, function(x) x[is_good_vert, ])  ## discard vertices with no BOLD
        }

        ## get training data and labels:

        B_train <- B[grep(train_session, names(B))]
        names(B_train) <- c("run1", "run2")
        ## initial downsampling -- discard trials that cannot be resampled :(
        B_train <- Map(
            function(data, discard) data[, which(!colnames(data) %in% discard)],
            data = B_train, discard = items_drop
        )
        ## split by stimulus set (pc50, bias):
        B_train <- list(
            pc50 = lapply(B_train, function(x) x[, colnames(x) %in% stimuli$pc50]),
            bias = lapply(B_train, function(x) x[, colnames(x) %in% stimuli$bias])
        )

        ## resample trial indices for training models:
        if (roi_i == 1) {  ## same order for all ROIs:        
            ## get trial indices:    
            trial_labels_pc50 <- lapply(B_train$pc50, colnames)
            trial_labels_bias <- lapply(B_train$bias, colnames)
            set.seed(0)
            downsamp_idx <- list(
                pc50 = lapply(trial_labels_pc50, resample_idx, n_resamples = n_resamples, resample_to = 3),
                bias = lapply(trial_labels_bias, resample_idx, n_resamples = n_resamples, resample_to = 3)
            )
            # ## build training response vectors:
            # y_train <- enlist(names(stimuli))
            # for (stimset in names(stimuli)) {
            #     train_labels <- c(colnames(downsamp_idx[[stimset]]$run1), colnames(downsamp_idx[[stimset]]$run2))
            #     y_train[[stimset]] <- cbind(
            #         td = train_labels,
            #         t = get_lower(train_labels),
            #         d = get_upper(train_labels),
            #         c = ifelse(tolower(get_upper(train_labels)) == get_lower(train_labels), "congr", "incon")
            #     )
            # }
        }

        
        res_stimset <- enlist(names(stimuli))
        for (stimset_i in seq_along(stimuli)) {

            stimset <- names(stimuli)[stimset_i]
            
            ## train model (baseline data) ----
            
            B_train_i <- B_train[[stimset]]
            B_train_c <- center_mom(B_train_i)  ## center at mean of means
            
            ## fit model:
            fits_t <- vector("list", n_resamples)
            fits_d <- vector("list", n_resamples)
            for (i in seq_len(n_resamples)) {

                training_set_i <- cbind(
                    B_train_c$run1[, downsamp_idx[[stimset]]$run1[i, ]], 
                    B_train_c$run2[, downsamp_idx[[stimset]]$run2[i, ]]
                    )
                #colnames(training_set_i) <- y_train[[stimset]][, "td"]

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
                    #preds_t <- lapply(fits_t, predict, newdata = t(scale2unit(B_test_t)), type = "posterior")
                    #preds_d <- lapply(fits_d, predict, newdata = t(scale2unit(B_test_d)), type = "posterior")
                    preds_t <- lapply(fits_t, mahals, newdata = t(scale2unit(B_test_t)))
                    preds_d <- lapply(fits_d, mahals, newdata = t(scale2unit(B_test_d)))
                } else if (decoder == "svm") {
                    preds_t <- lapply(fits_t, function(x) attr(predict(x, newdata = t(scale2unit(B_test_t)), probability = TRUE), "probabilities"))
                    preds_d <- lapply(fits_d, function(x) attr(predict(x, newdata = t(scale2unit(B_test_d)), probability = TRUE), "probabilities"))
                }
                
                ## summarize predictions

                ## aggregate predictions over resampled models:
                # preds_bar <- list(
                #     target = Reduce("+", lapply(preds_t, qlogis)) / n_resamples,
                #     distractor = Reduce("+", lapply(preds_d, qlogis)) / n_resamples
                # )                
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
                    "logit_target", "logit_distractor", "pred_target", "pred_distractor", 
                    "true_target", "true_distractor"
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


# a <- res
# lst <- list()
# for (i in seq(1, 191, 10)) {
#     lst <- c(lst, list(as.data.table(a[i:(i+9)])))
# }
# length(lst)
d <- rbindlist(res)

#saveRDS(d, here("out", "pda_lambda=1_pseudo=1_nresamp=1000_stimset-sep.RDS"))
#saveRDS(d, here("out", "pda_lambda=1_pseudo=2_nresamp=100_stimset-sep.RDS"))
saveRDS(d, here("out", "pda_lambda=100_pseudo=2_nresamp=100_stimset-sep_meas-mah.RDS"))

d_sum <- d %>%
    group_by(subj, roi, stimset) %>%
    summarize(
        cr_target = mean(pred_target == true_target),
        cr_distractor = mean(pred_distractor == true_distractor),
        logit_distractor = mean(-logit_distractor),
        logit_target = mean(-logit_target)
        )

d_suml <- d_sum %>%
    as.data.table %>%
    melt(
        id.vars = c("subj", "roi", "stimset"), 
        measure.vars = patterns(cr = "^cr", logit = "^logit"),
        variable.name = "model"
        )
d_suml$model <- c("target", "distractor")[d_suml$model]

d_suml %>%
    group_by(roi, stimset, model) %>%
    summarize(
        logit_t = t.test(logit)$statistic, cr_t = t.test(cr - 0.25)$statistic,
        cr = mean(cr), logit = mean(logit)
    ) %>%
    #filter(roi %in% fprois) %>%
    arrange(-logit_t) %>% View

d_sumsum <- d_suml %>%
    group_by(roi, stimset, model) %>%
    summarize(
        logit_p = t.test(logit, alternative = "greater")$p.value, cr_p = t.test(cr - 0.25, alternative = "greater")$p.value,
        logit_t = t.test(logit)$statistic, cr_t = t.test(cr - 0.25)$statistic,
        cr = mean(cr), logit = mean(logit)
    ) %>%
    group_by(stimset, model) %>%
    mutate(logit_p_fdr = p.adjust(logit_p, "fdr"), cr_p_fdr = p.adjust(cr_p, "fdr"))

d_sumsum %>%
    filter(logit_p_fdr < 0.05) %>% View

d_sumsum %>%
    filter(cr_p_fdr < 0.05) %>% View






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



