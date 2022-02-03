glmname <- "glmsingle_wave1"
subjlist <- "wave1_unrel"
sessions <- c("baseline", "proactive", "reactive")
waves <- "wave1"
tasks <- "Stroop"
subjs_rm <- ""#c("DMCC5820265", "DMCC9478705", "DMCC4260551")


library(here)
library(dplyr)
library(data.table)
library(gifti)
library(abind)
library(mfutils)
library(rhdf5)
library(doParallel)
library(foreach)
source(here("src", "stroop-rsa-pc.R"))
subjs <- fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]]
subjs <- setdiff(subjs, subjs_rm)
beh <- fread(here("in", "behavior-and-events_stroop_2021-10-20_nice.csv"))
beh <- beh[wave == "wave1" & subj %in% subjs]
beh$pc[beh$pc %in% c("mi", "mc")] <- "bias"
beh[, ttpc := paste0(trial_type, "_", pc)]
cols <- c("subj", "session", "run", "trial_num", "color", "word", "item", "trial_type", "pc", "ttpc", "time_target_onset")
beh <- beh[, ..cols]
## columns in beh that define the grouping variables for each decoding model:
varcols <- c(target = "color", distractor = "word", congruency = "ttpc")


## script-specific functions

onset_indicator <- function(onsets, .n_tr, .sec_tr = 1.2) {

  n_trial <- length(onsets)
  trgrid <- numeric(.n_tr)
  onsets_shifted <- onsets - onsets %% .sec_tr  ## assume that onsets are yoked with TR; moves back to begin w/ TR
  onsets_idx <- onsets_shifted / .sec_tr + 1  ## gives image (index for tr grid) to which onset belongs
  ## scaling by duration of the TR gives the index (which TR) it is, rel to 0 (add 1 b/c 1-based index needed)

  imat <- matrix(0, nrow = .n_tr, ncol = n_trial)
  imat[cbind(onsets_idx, seq_along(onsets_idx))] <- 1  ## mark onset of each trial (TR==1)
  
  imat

}

base_dir <- here("out", "glms")

## run ----


## iterators for dev:
session_i <- 1
subj_i <- 1
run_i <- 1


for (subj_i in seq_along(subjs)) {

    subj_nm <- subjs[subj_i]
    path_glm <- file.path(base_dir, subj_nm, "wave1", "RESULTS", "Stroop", glmname)
    if (!dir.exists(path_glm)) dir.create(path_glm)

    for (session_i in seq_along(sessions)) {

        session_nm <- sessions[session_i]
        X <- enlist(runs)

        for (run_i in 1:2) {
            
            run_nm <- runs[run_i]
            beh_i <- beh[session == session_nm & subj == subj_nm & run == run_i]
            beh_i <- beh_i[sort(trial_num), ]  ## just to make sure
            stopifnot(nrow(beh_i) == n_trial[session_nm])

            imat_onsets <- onset_indicator(onsets = beh_i$time_target_onset, .n_tr = n_tr[session_nm])
            imat_onsets1 <- onset_indicator1(onsets = beh_i$time_target_onset, .n_tr = n_tr[session_nm])
            Xi <- enlist(names(varcols))
            for (i in seq_along(Xi)) Xi[[i]] <- imat_onsets %*% indicator(beh_i[[varcols[i]]])

            X[[run_nm]] <- Xi

        }

        for (i in seq_along(varcols)) {
            
            Xi <- lapply(X, "[[", i)  ## all sessions and runs for given varcol
            stopifnot(length(unique(lapply(Xi, dim))) == 1)  ## dims must match across runs
            Xi <- abind(Xi, rev.along = 0)      
            file_name <- file.path(path_glm, paste0("design_", names(varcols[i]), ".h5"))
            h5write(Xi, file_name, session_nm)

        }

    }

}
