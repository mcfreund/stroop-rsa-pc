library(colorout)
library(here)
library(data.table)
source(here("src", "stroop-rsa-pc.R"))

s <- fread(here("out", "subjlist.csv"))

## development list: reactive test--retest
ispc_retest <- unique(s[is_ispc_retest == TRUE, "subj"])
fwrite(ispc_retest, here("out", "subjlist_ispc_retest.txt"), col.names = FALSE)
