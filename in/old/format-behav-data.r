library(colorout)
library(here)
library(data.table)

d <- rbindlist(
  list(
    wave1 = fread(here("in", "old", "dmcc2_behavior-and-events_stroop_2022-05-27.csv")),
    wave2 = fread(here("in", "old", "dmcc3_behavior-and-events_stroop_2022-05-27.csv")),
    wave3 = fread(here("in", "old", "dmcc4_behavior-and-events_stroop_2021-10-20.csv"))
  ),
  idcol = "wave"
)

d[session == "bas"]$session <- "baseline"
d[session == "rea"]$session <- "reactive"
d[session == "pro"]$session <- "proactive"

d$task <- "Stroop"

d$trial.type <- ifelse(d$trial.type == "i", "incon", "congr")

names(d) <- gsub("\\.", "_", names(d))

fwrite(d, here("in", "behavior-and-events_stroop_2022-05-27_nice.csv"))
