## about ----
##
## splits and moves movregs for run-wise GLMs.
##
## mike freund, 2020-04-04
## adapted for ub55 2020-08-20

library(colorout)
library(here)
library(dplyr)
library(magrittr)
library(data.table)

source(here("src", "evts", "write-events_funs.R"))

dir.analysis <- here("out", "glms")

s <- fread(here("out", "subjlist.csv"))
s <- s[needs_rerun == FALSE]  ## remove NA rows (why do these exist?)

to.split <- as.data.table(table(s$subj, s$session, s$wave))[N == 1]
names(to.split)[1:3] <- c("subj", "session", "wave")
to.split$task <- "Stroop"

to.split

movregs <- split.movregs(to.split, dir.to = dir.analysis)
movregs <- as.data.table(movregs)
sum(movregs$has.unexpected.nrow)
sum(movregs$is.missing.dir)
movregs[has.unexpected.nrow == TRUE]
movregs[is.missing.dir == TRUE]

fwrite(movregs, here("out", "glms", paste0("summary_split_movregs.csv")))
