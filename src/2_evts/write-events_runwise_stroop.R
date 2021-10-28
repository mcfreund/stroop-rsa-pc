## about ----
##
## writes events for run-wise GLMs.

library(colorout)
library(here)
library(dplyr)
library(magrittr)
library(data.table)
library(mikeutils)

source(here("src", "evts", "write-events_funs.R"))

d <- fread(here("in", "behavior-and-events_stroop_2021-10-20_nice.csv"))
d <- d[wave %in% c("wave1", "wave2")]

dir.analysis <- here("out", "glms")

s <- fread(here("out", "subjlist.csv"))



## run ----


## prep dfs for loops

d %<>% arrange(subj, session, trial_type)
d$time_target_onset_shift600 <- d$time_target_onset - 0.6  ## shift onset (s)
d$alltrials <- "alltrials"  ## make column for alltrials
l <- d %>% split(droplevels(interaction(.$subj, .$session, .$wave)))


## events

args_stroop_events <- expand.grid(
  var.level    = unique(d$item),
  name.var     = "item",
  dir.analysis = dir.analysis,
  name.onset   = "time_target_onset_shift600",
  by.run       = TRUE,
  fname.suffix = "_shift600",
  stringsAsFactors = FALSE
)

args_stroop_events

results_stroop_events <- lapply(l, write.events, .args = args_stroop_events)
results_stroop_events <- bind_rows(results_stroop_events, .id = "subj_session")
 
fwrite(results_stroop_events, here("out", "glms", paste0("summary_write-events_1rpm_shift600.csv")))
#unique(results_stroop_events$n.events)


## alltrials

args_stroop_alltrials <- expand.grid(
  var.level    = "alltrials",
  name.var     = "alltrials",
  dir.analysis = dir.analysis,
  name.onset   = "time_target_onset_shift600",
  by.run       = TRUE,
  fname.suffix = "_shift600",
  stringsAsFactors = FALSE
)

args_stroop_alltrials

results_stroop_alltrials <- lapply(l, write.events, .args = args_stroop_alltrials)
results_stroop_alltrials <- bind_rows(results_stroop_alltrials, .id = "subj_session")
fwrite(results_stroop_alltrials, here("out", "glms", paste0("summary_write-alltrials_1rpm_shift600.csv")))
