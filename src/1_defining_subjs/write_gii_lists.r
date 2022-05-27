library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

## first, source write_gii_lists.sh

fnames_subjlist <- list.files(here("in/subjs"), pattern = "^giis_dmcc*", full.names = TRUE)
subjs_fmri <- rbindlist(lapply(fnames_subjlist, fread, header = FALSE), idcol = "dmcc")
names(subjs_fmri)[2] <- "file"
subjs_fmri1 <- separate(subjs_fmri[dmcc == 1], file, c("dot", "subj", "dir1", "task", "session", "file"), sep = "/")
subjs_fmri2 <- separate(subjs_fmri[dmcc %in% 2:3], file, c("dir2", "dot", "subj", "dir1", "task", "session", "file"), sep = "/")

subjs_fmri <- rbind(
    subjs_fmri1[, c("dmcc", "subj", "session", "file")],
    subjs_fmri2[, c("dmcc", "subj", "session", "file")]
)
subjs_fmri$dmcc <- subjs_fmri$dmcc + 1
subjs_fmri$run <- ifelse(grepl("1_AP", subjs_fmri$file), "run1", ifelse(grepl("2_PA", subjs_fmri$file), "run2", NA))
subjs_fmri_sum <- subjs_fmri %>%
    group_by(dmcc, subj, session) %>%
    summarize(count = n()) %>%
    filter(count == 4) %>%
    select(dmcc, subj, session) %>%
    filter(!grepl("_", subj))

image(table(subjs_fmri_sum$subj, subjs_fmri_sum$session))

fwrite(subjs_fmri_sum, here("in", "subjs", paste0("subjlist_fmri_all_", Sys.Date(), ".txt")))