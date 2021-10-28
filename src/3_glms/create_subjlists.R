library(here)
s <- read.csv(here("out", "subjlist.csv"))
write.table(unique(s$subj), here("out", "glms", "allsubjs.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
