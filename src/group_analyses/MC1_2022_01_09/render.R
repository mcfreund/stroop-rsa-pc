library(here)
library(rmarkdown)
library(mfutils)

src_dir <- "group_analyses/MC1_2022_01_09"
roisets <- c("glasser2016_parcel", "schaefer2018_17_200_parcel")
ttype_subsets <- c("all")
measures <- c("crcor")
prewhs <- "none"

ps <- expand.grid(
  ttype_subset = ttype_subsets, prewh = prewhs, measure = measures, roiset = roisets,
  stringsAsFactors = FALSE
  )
for (i in seq_len(nrow(ps))) {
  p <- c(ps[i, ])
  render_report(name = "parcel", src_dir = src_dir, params = p)
}
