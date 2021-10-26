## about ----
##
## writes events for run-wise GLMs.
##
## mike freund, 2020-04-04
## adapted for ub55 2020-08-20

library(colorout)
library(here)
library(dplyr)
library(magrittr)
library(data.table)
library(mikeutils)

source(here("src", "evts", "write-events_funs.R"))

d <- rbindlist(
  list(
    wave1 = fread(here("in", "dmcc2_behavior-and-events_stroop_2021-10-20.csv")),
    wave2 = fread(here("in", "dmcc3_behavior-and-events_stroop_2021-10-20.csv"))
  ),
  idcol = "wave"
)

d[session == "bas"]$session <- "baseline"
d[session == "rea"]$session <- "reactive"
d[session == "pro"]$session <- "proactive"
d$task <- "Stroop"
d$trial.type <- ifelse(d$trial.type == "i", "incon", "congr")


dir.analysis <- here("out", "glms")

s <- fread(here("out", "subjlist.csv"))
s <- s[needs_rerun == FALSE]  ## remove NA rows (why do these exist?)


## format ----


## prep dfs for loops

d %<>% arrange(subj, session, trial.type)


## shift onset timing

d$time.target.onset.shift600 <- d$time.target.onset - 0.6  ## s


## split

l <- d %>% split(droplevels(interaction(.$subj, .$session, .$wave)))


## events

args.stroop.events <- expand.grid(
  var.level    = unique(d$item),
  name.var     = "item",
  dir.analysis = dir.analysis,
  name.onset   = "time.target.onset.shift600",
  by.run       = TRUE,
  fname.suffix = "_shift600",
  stringsAsFactors = FALSE
)

args.stroop.events

results.stroop.events <- lapply(l, write.events, .args = args.stroop.events)
results.stroop.events <- bind_rows(results.stroop.events, .id = "subj.session")
 

## write records

fwrite(results.stroop.events, here("out", "glms", paste0("summary_write-events_runwise_shift600.csv")))
#unique(results.stroop.events$n.events)