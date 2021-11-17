library(colorout)
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
