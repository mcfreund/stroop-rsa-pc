args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

glmname <- args$glmname
atlas_name <- args$atlas_name
roi_col <- args$roi_col
subjlist <- args$subjlist
waves <- strsplit(args$waves, "\\|")[[1]]
sessions <- strsplit(args$sessions, "\\|")[[1]]
subjects <- data.table::fread(here("out", "subjs", paste0("subjlist_", subjlist, ".txt")))[[1]]
space <- args$space

if ("n_cores" %in% names(args)) n_cores <- as.numeric(args$n_cores)  ## optional arg that needs coersion
if ("n_resamples" %in% names(args)) n_resamples <- as.numeric(args$n_resamples)
if ("prewh" %in% names(args)) prewh <- args$prewh
if ("measure" %in% names(args)) measure <- args$measure
if ("ttype_subset" %in% names(args)) ttype_subset <- args$ttype_subset
if ("expected_min" %in% names(args)) expected_min <- as.numeric(args$expected_min)
if ("overwrite" %in% names(args)) overwrite <- args$overwrite == "TRUE"
if ("delete_files" %in% names(args)) delete_files <- args$delete_files == "TRUE"
if ("suffix" %in% names(args)) suffix <- args$suffix