## input arguments

n_cores <- 4
atlas_roi_nms <- c("glasser2016__superparcel", "glasser2016__parcel", "schaefer2018_17_200__parcel")
suffix <- "__nresamp100__lsall_1rpm__divnormed__demeaned"
waves <- c("wave1", "wave2")

## libraries

library(here)
library(dplyr)
library(tidyr)
library(mfutils)
library(data.table)
library(foreach)
library(doParallel)

source(here("src", "stroop-rsa-pc.R"))

## constants

subset_nms <- c("bias", "pc50", "biaspc50")
names(ttypes)[names(ttypes) == "all"] <- "biaspc50"
path_data <- here("out", "decoding")
subjdirs <- list.dirs(path_data)[-1]


## run

for (wave in waves) {

  for (atlas_roi_nm in atlas_roi_nms) {

    ## filename information:

    fname_in <- paste0("trainbaseline__", wave, "__svm__", atlas_roi_nm, suffix, ".RDS")
    fname_out <- paste0("trainbaseline__", wave, "__svm__", atlas_roi_nm, suffix, "_summarized.RDS")

    does_exist <- file.exists(file.path(subjdirs, fname_in))
    subjdirs <- subjdirs[does_exist]

    ## loop over input files

    cl <- makeCluster(n_cores, type = "FORK")
    registerDoParallel(cl)
    foreach(subjdir_i = seq_along(subjdirs), .inorder = FALSE, .verbose = TRUE) %dopar% {

      subjdir <- subjdirs[subjdir_i]

      d <- readRDS(file.path(subjdir, fname_in))

      ## add, remove relevant cols:
      d[, 1 := NULL]  ## remove duplicate column "roi"
      d[, subj := gsub(".*decoding/", "", subjdir)]
      d[, target := get_variable(stimulus, "target")]
      d[, distractor := get_variable(stimulus, "distractor")]
      d[, congruency := get_variable(stimulus, "congruency")]

      ## extract looping variables
      nms <- names(d)
      rois <- unique(d$roi)
      sessions <- unique(d$session)

      ## loop over different decoders and calculate "evidence"
      id <- combo_paste(subset_nms, c("target", "distractor", "congruency"), rois, sessions, sep = "__")
      results <- enlist(id)
      for (ttype_subset in subset_nms) {
        di <- d[train_stimset == ttype_subset]
        for (var_name in c("target", "distractor", "congruency")) {  ## loop over factors
          if (var_name != "congruency" & ttype_subset == "biaspc50") next
          dii <- di[variable == var_name]
          for (roi_name in rois) {
            diii <- dii[roi == roi_name]
            for (session_name in sessions) {


              d_ses <- diii[session == session_name]
              ## this builds a list of "decoders": each list element corresponds to a single binary classification/
              ## distinction,
              ## consisting of a character vector with the class names within the distinction
              ## the first element of the character vector corresponds to the class assigned a positive sign
              all_classes <- unique(d_ses[[var_name]])
              var_cols <- nms[grepl(paste0("proj.*", all_classes, collapse = "|"), nms)]  ## get colnms of relevant prjs
              raw_projs <- d_ses[, ..var_cols]  ## get just cols with data
              var_cols <- var_cols[!is_equal(colMeans(is.na(raw_projs)), 1)]  ## take only colnms that have observations
              decoders <- strsplit(gsub("proj.", "", var_cols), "/")
              names(decoders) <- var_cols

              ## now loop over each binary decoder and transform decision values
              res_decoders <- enlist(names(decoders))
              for (decoder_i in seq_along(decoders)) {

                classes <- decoders[[decoder_i]]  ## pos, neg
                colnm <- names(decoders)[decoder_i]  ## get relevant colnames

                ## get relevant data:
                is_class <- d_ses[[var_name]] %in% classes  ## levels (particular distinction)
                d_ses_i <- d_ses[is_class, ]

                ## center at mean of class means:
                y <- rbind(d_ses_i[[colnm]])  ## extract data
                g <- d_ses_i[[var_name]]  ## class labels
                mubar <- mean(average(y, g))
                y_c <- c(y - mubar)

                ## multiply by sign assigned to each class (i.e., bring negative class positive):
                class_sign <- sign((g == classes[1]) - 0.5)
                evidence <- y_c * class_sign

                ## bind in datatable
                newcol <- paste0(var_name, "_", paste0(classes, collapse = "_"))
                res <- data.table(evidence)
                #names(res) <- newcol
                res[, trial := d_ses_i$trial]  ## add trial info for binding
                res[,
                 c("roi", "session", "variable", "train_stimset") := 
                  list(roi_name, session_name, var_name, ttype_subset)
                 ]

                ## return in list:
                res_decoders[[decoder_i]] <- res

              }

              ##  for each trial, average over projections onto different decoders relevant for that trial
              res_decoders_sum <- rbindlist(res_decoders, idcol = "decoder")
              res_decoders_sum <- res_decoders_sum[
                  , .(evidence = mean(evidence)),
                  by = c("trial", "roi", "session", "variable", "train_stimset")
                  ]

              ## return in list:
              iter_id <- paste(ttype_subset, var_name, roi_name, session_name, sep = "__")
              results[[iter_id]] <- res_decoders_sum


            }
          }
        }
      }

      ## prepare relevant cols of input dataframe for binding to output
      cols_keep <- names(d)[grep("proj", names(d))]
      d[, (cols_keep) :=  NULL]  ##  remove now-extraneous cols

      ## merge output data to relevant cols of input data
      r <- rbindlist(results)
      output <- merge(d, r)

      ## save
      saveRDS(output, file.path(subjdir, fname_out))
      print(file.path(subjdir, fname_out))

    }
    stopCluster(cl)

    print(paste0(wave, " ", atlas_roi_nm))

  }

}