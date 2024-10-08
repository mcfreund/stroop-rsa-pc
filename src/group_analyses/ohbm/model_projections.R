## libraries

library(here)
library(knitr)
library(ggplot2)
library(mfutils)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)
library(broom)
library(parallel)
library(dplyr)
library(tidyr)
library(data.table)
library(brms)

source(here("src", "stroop-rsa-pc.R"))

## input arguments

n_chains <- 4
n_cores_outer <- 5
wave <- "wave1"
atlas_roi_nm <- "glasser2016__parcel"
# rois <-
#   c(
#     mfutils::combo_paste(c("L_", "R_"), c("V1", "V2", "V3")),
#     glasser2016_fsaverage5$key %>% filter(network == "Language", grepl("L_", parcel)) %>% pull(parcel)
#   )
rois <- mfutils::combo_paste(c("L_", "R_"), md$core)
stimsets <- c("bias")
#stimsets <- c("bias", "pc50", "biaspc50")
#atlas_roi_nm <- "schaefer2018_17_200__parcel"
#atlas_roi_nms <- c("glasser2016__superparcel", "glasser2016__parcel", "schaefer2018_17_200__parcel")
suffix <- "__nresamp100__lsall_1rpm__divnormed__demeaned__centered"
rois <- c("R_SCEF", "R_IFJp", "R_8BM", "L_V3", "L_V2", "L_STV", "L_STSdp")

## check number of cores allowable:

stopifnot(n_cores_outer*n_chains <= detectCores())

## utilities

fit_models <- function(
  x, formula, n_cores_outer, n_chains,
  return_fixef = TRUE,
  dir_out = path_out,
  ...
  ) {

  dir.create(dir_out, showWarnings = FALSE)

  mclapply(
    names(x),
    function(
      nm, 
      .data = x[[nm]],
      .formula = formula,
      .chains = n_chains,
      .return_fixef = return_fixef,
      .dir_out = dir_out,
      .file_name = fname,
       ...
      ) {

      tryCatch(
        {
          fit <- brm(formula = .formula, data = .data, chains = .chains, cores = .chains, ...)
          saveRDS(fit, file.path(dir_out, nm))
          if (.return_fixef) fixef(fit)
        },
        error = NA
      )

    },
    mc.cores = n_cores_outer,
    ...
  )
}


## constants

models <- c("target", "distractor", "congruency")
path_data <- here("out", "decoding")
subjdirs <- list.dirs(path_data)[-1]
path_data <- here("out", "decoding")
file_name <- paste0("trainbaseline__", wave, "__svm__", atlas_roi_nm, suffix, "_summarized")
does_exist <- file.exists(file.path(subjdirs, paste0(file_name, ".RDS")))
path_out <- here("out", "decoding", "_group", file_name)
subjdirs <- subjdirs[does_exist]
subjs <- gsub(".*decoding/", "", subjdirs)

subset_nms <- c("bias", "pc50", "biaspc50")
names(ttypes)[names(ttypes) == "all"] <- "biaspc50"



## data

l <- mclapply(subjdirs, function(x) readRDS(file.path(x, paste0(file_name, ".RDS"))), mc.cores = n_cores_outer*n_chains)
d <- rbindlist(l)
d <- d[roi %in% rois]

## add behavior:
beh <- fread(here("in", "behavior-and-events_stroop_2022-05-27_nice.csv"))
beh <- beh %>%
  group_by(subj, session, wave, run) %>%
  mutate(trial = trial_num + ifelse(session %in% c("baseline", "proactive"), 108, 120) * (run - 1)) %>%
  select(wave, session, run, trial, subj, time_target_onset, rt, acc)
dat <- data.table::merge.data.table(d, beh, by = c("subj", "session", "wave", "trial"), all.x = TRUE)
stopifnot(nrow(dat) == nrow(d))

## reformat cols:
dat[, region := gsub("L_|R_", "", roi)]
dat[, hemi := ifelse(grepl("L_", roi), "left", "right")]
dat[, region_hemi := paste0(toupper(substr(hemi, 1, 1)), "_", region)]
dat <- dat[!(region %in% c("10pp") & hemi == "left")]  ## not in glasser ggseg atlas?
## add test_stimset:
dat[, test_stimset := ""]
dat[stimulus %in% ttypes$bias, test_stimset := "bias"]
dat[stimulus %in% ttypes$pc50, test_stimset := "pc50"]

## code regressors

dat[, run := paste0("run", run)]


## TODO: subject lists, twins, etc.


## brms models

f_session <- brmsformula(
  evidence ~ 0 + session + run + (0 + session | subj),
  sigma ~ 0 + session + run + (0 + session | subj)
  )

f_session_congruency <- brmsformula(
  evidence ~ 0 + interaction(session, congruency) + (0 + interaction(session, congruency) | subj),
  sigma ~ 0 + interaction(session, congruency) + (0 + interaction(session, congruency) | subj)
  )

f_session_rtarget <- brmsformula(
  evidence ~ 0 + session + (0 + session | subj) + (0 + session | target),
  sigma ~ 0 + session + (0 + session | subj)
  )



## fit session model


for (model in models) {
  for (train_stimset_nm in stimsets) {

    if (train_stimset_nm == "biaspc50") {
      if (model != "congruency") next
      test_stimset_nm <- "pc50"
    } else {
      test_stimset_nm <- train_stimset_nm
    }

    print(c(model, train_stimset_nm, test_stimset_nm))

    l <-
      dat[variable == model & train_stimset == train_stimset_nm & test_stimset == test_stimset_nm] %>%
      split(.$region_hemi)

    fit_models(
      l,
      formula = f_session,
      n_cores_outer = n_cores_outer,
      n_chains = n_chains,
      family = student,
      dir_out = file.path(path_out, paste0(model, "_", train_stimset_nm))
      )

  }
}
