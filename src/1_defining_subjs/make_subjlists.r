library(here)
library(data.table)
source(here("src", "stroop-rsa-pc.R"))

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
