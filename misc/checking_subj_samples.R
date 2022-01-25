library(here)
library(data.table)

subjs <- fread(here("out", "subjlist_mcmi.txt"))[[1]]
subjs_mc1 <- fread(here("out", "subjlist_mc1_unrel.txt"))[[1]]
subjs_mc11 <- fread(here("out", "subjlist_mc1_unrel1.txt"))[[1]]
#subjs_mc1 <- unique(dat$subject)

length(subjs)
length(subjs_mc1)
length(subjs_mc11)

length(setdiff(subjs, subjs_mc1))
length(setdiff(subjs_mc1, subjs))

setdiff(subjs_mc11, subjs_mc1)

subjs_mz <- data.table(
  a = c(
    "DMCC1624043", "DMCC1971064", "233326", "209228", "DMCC4260551", "DMCC6755891", "DMCC6904377", "DMCC9850294",
    "162026", "DMCC8078683", "DMCC7297690", "DMCC6960387", "DMCC5065053", "DMCC5820265", "DMCC3206338"
  ),
  b = c(
    "DMCC3204738", "DMCC3509558", "352738", "765864", "DMCC7921988", "DMCC8050964", "DMCC9953810", "DMCC2560452",
    "568963", "DMCC6484785", "DMCC1596165", "DMCC8260571", "DMCC4368773", "DMCC3091953", "DMCC4854984"
  )
)
subjs_mz_mbsr <- data.table(
  a = c(
    "182436", "115825", "250427", "171330", "DMCC2834766", "DMCC2803654", "DMCC4191255", "562345",
    "DMCC5195268", "DMCC6661074", "DMCC5775387", "DMCC6371570", "DMCC2442951", "DMCC1328342",
    "DMCC3963378", "DMCC3876181", "DMCC5244053"
  ),
  b = c(
    "178647", "178243", "877168", "393550", "DMCC6671683", "DMCC9478705", "DMCC6418065", "130518",
    "DMCC8033964", "DMCC6705371", "DMCC6721369", "DMCC9441378", "DMCC6627478", "DMCC5009144",
    "DMCC2609759", "DMCC0472647", "DMCC8760894"
  )
  
)
subjs_dz <- data.table(
  a = c("198855", "DMCC8214059"),
  b = c("623844", "DMCC3062542")
)
subjs_dz_mbsr <- data.table(
  a = c("130114"),
  b = c("155938")
)
subjs_twins <- rbindlist(list(mz = subjs_mz, mz = subjs_mz_mbsr, dz = subjs_dz, dz = subjs_dz_mbsr), idcol = "twin")
subjs_twins[, twinpair := 1:.N]
subjs_twins <- melt(subjs_twins, id.vars = c("twin", "twinpair"), value.name = "subj")[, -"variable"]


setdiff(setdiff(subjs, subjs_mc1), subjs_twins$subj)
setdiff(setdiff(subjs_mc1, subjs), subjs_twins$subj)


