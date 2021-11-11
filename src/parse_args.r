args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

glmname <- args$glmname
roiset <- args$roiset
subjlist <- args$subjlist
waves <- strsplit(args$waves, "\\|")[[1]]
sessions <- strsplit(args$sessions, "\\|")[[1]]
if ("n_cores" %in% names(args)) args$n_cores <- as.numeric(args$n_cores)  ## optional arg that needs coersion

subjects <- data.table::fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]]

lapply(names(args), function(x, l) assign(x, l[[x]]), l = args)  ## assign to variable names in environment

## check some vals
if (exists("prewh")) stopifnot(prewh %in% c("obsall", "obsresamp"))
if (exists("measure")) stopifnot(measure %in% c("cveuc", "crcor"))
