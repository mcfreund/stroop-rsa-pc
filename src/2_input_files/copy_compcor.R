## about ----
##
## splits and moves movregs for run-wise GLMs.
##
## mike freund, 2020-04-04
## copies fmriprep compcor timeseries for use in afni GLMs

library(here)
library(dplyr)
library(data.table)


dir_analysis <- here("out", "glms")

s <- fread(here("out", "subjlist.csv"))

input <- as.data.table(table(s$subj, s$session, s$wave))[N == 1]
names(input)[1:3] <- c("subj", "session", "wave")
input$task <- "Stroop"
input$ses <- substr(input$session, 1, 3)
input$run1 <- "AP_run-1"
input$run2 <- "PA_run-2"
input <- melt(input, value.name = "run", measure.vars = c("run1", "run2"))
input$variable <- NULL


## build paths
dir_fmriprep_preprocessed <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_PREPROCESSED/"
input$dir_confounds <- file.path(
    dir_fmriprep_preprocessed,
    input$subj, 
    "derivatives", 
    "fmriprep", 
    paste0("sub-", input$subj), 
    paste0("ses-", input$wave, input$ses), 
    "func"
)
input$file_name <- 
    paste0(
        "sub-", input$subj, "_ses-", input$wave, input$ses, "_task-Stroop_acq-mb4", input$run, 
        "_desc-confounds_regressors.tsv"
        )

input$dir_exists <- sapply(input$dir_confounds, dir.exists)
input$file_exists <- file.exists(file.path(input$dir_confounds, input$file_name))
input$full_path <- file.path(input$dir_confounds, input$file_name)

input_existing <- input[dir_exists == TRUE & file_exists == TRUE]

table(input_existing$subj, paste0(input_existing$session, input_existing$wave))  ## only wave 1?

res <- parallel::mclapply(
    input_existing$full_path,
    function(x) {
        A <- fread(x)
        nms <- grep("a_comp_cor", names(A))
        A[, ..nms]
    },
    mc.cores = 12
)
unique(sapply(res, ncol))  ## all have 6 columns.
unique(sapply(res, nrow))  ## all have 540 or 580 nrows (+/-1)

## write

input_existing$run_rename <- gsub("AP_run-|PA_run-", "run", input_existing$run)
input_existing$out_file_name <- here(
    "out", "glms", input_existing$subj, "wave1", "INPUT_DATA", "Stroop", input_existing$session, 
    paste0("a_comp_cor_00-06_", input_existing$run_rename, ".1D")
    )
Map(
    function(x, y) fwrite(x, y, sep = " ", col.names = FALSE),
    x = res,
    y = input_existing$out_file_name
)
