args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

glmname <- args$glmname
roiset <- args$roiset
subjlist <- args$subjlist
waves <- strsplit(args$waves, "\\|")[[1]]
sessions <- strsplit(args$sessions, "\\|")[[1]]
subjects <- data.table::fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]]

if ("n_cores" %in% names(args)) n_cores <- as.numeric(args$n_cores)  ## optional arg that needs coersion
if ("prewh" %in% names(args)) prewh <- args$prewh
if ("measure" %in% names(args)) measure <- args$measure
if ("ttype_subset" %in% names(args)) ttype_subset <- args$ttype_subset
if ("expected_min" %in% names(args)) expected_min <- as.numeric(args$expected_min)
if ("overwrite" %in% names(args)) overwrite <- args$overwrite == "TRUE"