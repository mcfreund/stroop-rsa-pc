library(here)
library(rmarkdown)
library(mfutils)
library(reticulate)

src_dir <- "group_analyses/MCvMI_parcel200_2021-12-20"
ttype_subsets <- c("bias", "pc50")
prewhs <- "none"
measures <- "crcor"

ps <- expand.grid(ttype_subset = ttype_subsets, prewh = prewhs, measure = measures, stringsAsFactors = FALSE)
for (i in seq_len(nrow(ps))) {
  p <- c(ps[i, ])
  render_report(name = "parcel200", src_dir = src_dir, params = p)
}
