library(here)
library(data.table)
library(dplyr)

source(here("src", "stroop-rsa-pc.R"))
sample_cotwins <- function(x) {
  set.seed(0)
  twinpairs <- unique(x$twinpair[duplicated(x$twinpair) & !is.na(x$twinpair)])
  x <- x %>% filter(twinpair %in% twinpairs)  ## only twins in sample
  x %>%
    group_by(twinpair) %>%
    sample_n(1) %>%
    pull(subj)
}

snew <- fread(here("out", paste0("subjlist_2022-05-27.csv")))
## single session lists
glmname <- "lsall_1rpm"
for (wav in c("wave1", "wave2")) {
  for (ses in c("baseline", "proactive", "reactive")) {
    subjs <- snew[session == ses & wave == wav, subj]
    fnames <- here(
      "out", "glms", subjs, wav, "RESULTS", "Stroop", 
      paste0(ses, "_", glmname), 
      paste0("STATS_", subjs, "_run1_L_REML.func.gii")
        )
    file_exists <- file.exists(fnames)
    print(paste0(sum(file_exists), "/", length(file_exists)))
    print(c("missing:", subjs[!file_exists]))
    subjs <- as.data.frame(subjs[file_exists])
    fname <- here("out", "subjs", paste0("subjlist_", ses, "_", wav, ".txt"))
    fwrite(subjs, fname, col.names = FALSE)
  }
}
## cross-session lists
for (wav in c("wave1", "wave2")) {
  
    subjs <- snew[wave == wav, unique(subj)]
    #is_good_subj <- table(subjs) == 3
    subjs_fitted <- c()
    for (ses in sessions) {
      fnames <- here(
        "out", "glms", subjs, wav, "RESULTS", "Stroop", 
        paste0(ses, "_", glmname), 
        paste0("STATS_", subjs, "_run1_L_REML.func.gii")
          )
      file_exists <- file.exists(fnames)
      subjs_fitted <- c(subjs_fitted, subjs[file_exists])
    }
    counts <- table(subjs_fitted)
    subjs_complete <- names(counts)[counts == 3]
    print(length(subjs_complete))    
    
    fname <- here("out", "subjs", paste0("subjlist_", wav, ".txt"))
    fwrite(as.data.frame(subjs_complete), fname, col.names = FALSE)

}






s <- fread(here("out", "subjlist.csv"))

## wave 1 list: baseline, proactive, reactive

subjs_wave1 <- names(which(table(s[wave == "wave1", subj]) == 3))
s_wave1 <- s[wave == "wave1" & subj %in% subjs_wave1]
wave1_cotwins_exclude <- sample_cotwins(s_wave1)
s_wave1_unrel <- unique(s_wave1[!subj %in% wave1_cotwins_exclude, "subj"])
fwrite(s_wave1_unrel, here("out", "subjlist_wave1_unrel.txt"), col.names = FALSE)
fwrite(unique(s_wave1[, "subj"]), here("out", "subjlist_wave1.txt"), col.names = FALSE)


## development list: reactive test--retest

ispc_retest <- unique(s[is_ispc_retest == TRUE, "subj"])
fwrite(ispc_retest, here("out", "subjlist_ispc_retest.txt"), col.names = FALSE)


all_retest <- intersect(
    unique(s[is_ispc_retest == TRUE | is_ispc_retest_cotwin == TRUE, subj]), 
    unique(s[is_mc_mi_retest == TRUE | is_mc_mi_retest_cotwin == TRUE, subj])
    )
fwrite(as.data.table(all_retest), here("out", "subjlist_all_retest.txt"), col.names = FALSE)


## primary analysis list: mcmi (wave1) and mimc (wave1/wave2)

mcmi <- unique(s[is_mc_mi == TRUE, "subj"])
mimc <- unique(s[is_mi_mc == TRUE, "subj"])
mcmi_cotwin <- unique(s[is_mc_mi_cotwin == TRUE, "subj"])
mimc_cotwin <- unique(s[is_mi_mc_cotwin == TRUE, "subj"])

fwrite(mcmi, here("out", "subjlist_mcmi.txt"), col.names = FALSE)
fwrite(mimc, here("out", "subjlist_mimc.txt"), col.names = FALSE)
fwrite(mcmi_cotwin, here("out", "subjlist_mcmi_cotwin.txt"), col.names = FALSE)
fwrite(mimc_cotwin, here("out", "subjlist_mimc_cotwin.txt"), col.names = FALSE)

## per session*wave

fwrite(s[session == "baseline" & wave == "wave1", "subj"], here("out", "subjlist_mc1.txt"), col.names = FALSE)
fwrite(s[session == "baseline" & wave == "wave2", "subj"], here("out", "subjlist_mc2.txt"), col.names = FALSE)
fwrite(s[session == "proactive" & wave == "wave1", "subj"], here("out", "subjlist_mi1.txt"), col.names = FALSE)
fwrite(s[session == "proactive" & wave == "wave2", "subj"], here("out", "subjlist_mi2.txt"), col.names = FALSE)


## jneurosci (re)analysis list: wave1 proactive

strooprsa <- fread("https://raw.githubusercontent.com/mcfreund/stroop-rsa/master/in/behavior-and-events_group201902.csv")
strooprsa <- unique(strooprsa[, .(subj, is.analysis.group)])
strooprsa[, .N, by = is.analysis.group]
reanalysis_primary <- strooprsa[is.analysis.group == TRUE, "subj"]
reanalysis_cotwin <- strooprsa[is.analysis.group == FALSE, "subj"]

fwrite(reanalysis_primary, here("out", "subjlist_jneurosci_primary.txt"), col.names = FALSE)
fwrite(reanalysis_cotwin, here("out", "subjlist_jneurosci_cotwin.txt"), col.names = FALSE)

reanalysis <- c(reanalysis_primary$subj, reanalysis_cotwin$subj)

## jneurosci replication list: wave1 proactive NOT IN reanalysis

replication <- s[session == "proactive" & wave == "wave1" & !subj %in% reanalysis]
replication_cotwins <- sample_cotwins(replication)
replication_primary <- replication[!subj %in% replication_cotwins, subj]
length(replication_cotwins)
length(replication_primary)
# intersect(c(replication_primary, replication_cotwins), c(reanalysis_primary, reanalysis_cotwin))
# intersect(replication_primary, replication_cotwins)
# table(replication[subj %in% replication_primary, twinpair])
# table(replication[subj %in% replication_cotwins, twinpair])

fwrite(as.data.frame(replication_primary), here("out", "subjlist_jneurosci_replication_primary.txt"), col.names = FALSE)
fwrite(as.data.frame(replication_cotwins), here("out", "subjlist_jneurosci_replication_cotwin.txt"), col.names = FALSE)


## MC1 list
s_mc1 <- s[session == "baseline" & wave == "wave1"]
mc1_cotwins_exclude <- sample_cotwins(s_mc1)
fwrite(s_mc1[!subj %in% mc1_cotwins_exclude, "subj"], here("out", "subjlist_mc1_unrel1.txt"), col.names = FALSE)



