library(here)
library(rmarkdown)
library(mfutils)

src_dir <- "group_analyses/MC1_2022_01_09"
roisets <- c("glasser2016_parcel", "schaefer2018_17_200_parcel")
measures <- "crcor"
prewhs <- "none"

ps <- expand.grid(
  prewh = prewhs, measure = measures, roiset = roisets,
  stringsAsFactors = FALSE
  )
for (i in seq_len(nrow(ps))) {
  p <- c(ps[i, ])
  render_report(name = "parcel_pc50-bias-contrast", src_dir = src_dir, params = p)
}
