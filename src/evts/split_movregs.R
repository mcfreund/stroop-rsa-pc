## about ----
##
## splits and moves movregs for run-wise GLMs.
##
## mike freund, 2020-04-04
## adapted for ub55 2020-08-20

library(here)
library(dplyr)
library(magrittr)
library(data.table)

source(here("src", "evts", "write-events_funs.R"))

dir.analysis <- here("out", "glms")

s <- fread(here("out", "subjlist.csv"))


to.split <- expand.grid(
  subj = subjs,
  task = c("Stroop"),
  session = c("baseline", "proactive", "reactive"),
  wave = c("wave1", "wave2"),
  stringsAsFactors = FALSE
)

movregs <- split.movregs(to.split, dir.to = dir.analysis)
sum(movregs$has.unexpected.nrow)
sum(movregs$is.missing.dir)


