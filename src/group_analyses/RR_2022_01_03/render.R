library(here)
library(rmarkdown)
library(mfutils)

src_dir <- "group_analyses/RR_2022_01_03"
roisets <- c("glasser2016_parcel", "schaefer2018_17_200_parcel")
ttype_subsets <- c("bias", "pc50")
measures <- c("cveuc", "crcor")
#prewhs <- c("none", "obsall")
prewhs <- "obsallave"

ps <- expand.grid(
  ttype_subset = ttype_subsets, prewh = prewhs, measure = measures, roiset = roisets,
  stringsAsFactors = FALSE
  )
for (i in seq_len(nrow(ps))) {
  p <- c(ps[i, ])
  if (p$measure == "crcor" & p$prewh == "obsallave") next
  render_report(name = "parcel", src_dir = src_dir, params = p)
}
