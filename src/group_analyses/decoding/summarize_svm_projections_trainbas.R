## input arguments

wave <- "wave1"
atlas_roi_nms <- c("glasser2016__superparcel", "glasser2016__parcel", "schaefer2018_17_200__parcel")
suffix <- "__nresamp100__lsall_1rpm__divnormed__demeaned"

## libraries

library(here)
library(data.table)
library(knitr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mfutils)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)

source(here("src", "stroop-rsa-pc.R"))

## constants

subset_nms <- c("bias", "pc50", "biaspc50")
names(ttypes)[names(ttypes) == "all"] <- "biaspc50"
path_data <- here("out", "decoding")
subjdirs <- list.dirs(path_data)[-1]


## run

for (atlas_roi_nm in atlas_roi_nms) {

  ## filename information:

  fname_in <- paste0("trainbaseline__", wave, "__svm__", atlas_roi_nm, suffix, ".RDS")
  fname_out <- paste0("trainbaseline__", wave, "__svm__", atlas_roi_nm, suffix, "_summarized.RDS")

  does_exist <- file.exists(file.path(subjdirs, fname_in))
  subjdirs <- subjdirs[does_exist]

  ## loop over input files

  for (subjdir in subjdirs) {

    d <- readRDS(file.path(subjdir, fname_in))

    ## add relevant colnames
    d$target <- get_variable(d$stimulus, "target")
    d$distractor <- get_variable(d$stimulus, "distractor")
    d$congruency <- get_variable(d$stimulus, "congruency")
    d <- d[, -1]
    d$subj <- gsub(".*decoding/", "", subjdir)

    ## extract looping variables
    nms <- names(d)
    rois <- unique(d$roi)
    sessions <- unique(d$session)

    ## initialize results data frame:
    d_new <- d
    cols_keep <- names(d)[-grep("proj", names(d_new))]
    d_new <- d_new[, ..cols_keep]

    ## loop over different decoders and calculate "evidence"

    for (ttype_subset in subset_nms) {
      for (variable in c("target", "distractor", "congruency")) {  ## loop over factors

        if (variable != "congruency" & ttype_subset == "biaspc50") next

        for (roi in rois) {
          for (session in sessions) {

            ## this builds a list of "decoders": each list element corresponds to a single binary classification/
            ## distinction,
            ## consisting of a character vector with the class names within the distinction
            ## the first element of the character vector corresponds to the class assigned a positive sign
            all_classes <- unique(d[[variable]])
            var_cols <- nms[grepl(paste0("proj.*", all_classes, collapse = "|"), nms)]  ## get colnames of relevant prjs
            is_val_decoder <- 
              d$roi %in% roi &
              d$session %in% session &
              d$variable == variable &  ## factor
              d[, train_stimset] == ttype_subset  ## training/testing set
            raw_projs <- d[is_val_decoder, ..var_cols]  ## get just cols with data
            var_cols <- var_cols[!is_equal(colMeans(is.na(raw_projs)), 1)]  ## take only colnames that have observations
            decoders <- strsplit(gsub("proj.", "", var_cols), "/")
            names(decoders) <- var_cols

            ## now loop over each binary decoder and transform decision values

            for (decoder_i in seq_along(decoders)) {

              classes <- decoders[[decoder_i]]  ## pos, neg

              ## get relevant data:
              is_val <-
                d$roi %in% roi &
                d$session %in% session &
                d$variable == variable &  ## factor
                d[[variable]] %in% classes &  ## levels (particular distinction)
                d[, train_stimset] == ttype_subset  ## training/testing set
              colnm <- names(decoders)[decoder_i]  ## get relevant colnames
              di <- d[is_val, .SD, .SDcols = c(variable, colnm)]

              ## center at mean of class means:
              y <- rbind(di[[colnm]])  ## extract data
              g <- di[[variable]]  ## class labels
              mubar <- mean(average(y, g))
              y_c <- c(y - mubar)

              ## multiply by sign assigned to each class (i.e., bring negative class positive):
              class_sign <- sign((g == classes[1]) - 0.5)
              evidence <- y_c * class_sign

              ## return to datatable:
              newcol <- paste0(variable, "_", paste0(classes, collapse = "_"))
              d_new[is_val, (newcol) := evidence]

            }
          }
        }
      }
    }

    ## now summarize across different decoders (for multiclass classifiers)

    nms_target <- names(d_new)[grepl("target_", names(d_new))]
    nms_distractor <- names(d_new)[grepl("distractor_", names(d_new))]
    nms_congruency <- names(d_new)[grepl("congruency_", names(d_new))]
    d_new[, proj_target := rowMeans(.SD, na.rm = TRUE), .SDcols = nms_target]
    d_new[, proj_distractor := rowMeans(.SD, na.rm = TRUE), .SDcols = nms_distractor]
    d_new[, proj_congruency := rowMeans(.SD, na.rm = TRUE), .SDcols = nms_congruency]

    ## remove now-extraneous projection cols
    cols_rm <- c(nms_target, nms_distractor, nms_congruency)
    cols_keep_final <- names(d_new)[!names(d_new) %in% cols_rm]
    d_new <- d_new[, ..cols_keep_final]


    saveRDS(d_new, file.path(subjdir, fname_out))
    print(file.path(subjdir, fname_out))


  }
}