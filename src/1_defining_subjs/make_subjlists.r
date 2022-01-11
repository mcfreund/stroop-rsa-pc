library(here)
library(data.table)
library(dplyr)

source(here("src", "stroop-rsa-pc.R"))
sample_cotwins <- function(x) {
  set.seed(0)
  twinpairs <- unique(x$twinpair[duplicated(x$twinpair)])
  x <- x %>% filter(twinpair %in% twinpairs)  ## only twins in sample
  x %>%
    group_by(twinpair) %>%
    sample_n(1) %>%
    pull(subj)
}

s <- fread(here("out", "subjlist.csv"))

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
fwrite(s_mc1[!subj %in% mc1_cotwins_exclude, "subj"], here("out", "subjlist_mc1_unrel.txt"), col.names = FALSE)
