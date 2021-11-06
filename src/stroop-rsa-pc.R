## settings ----

options(datatable.print.trunc.cols = TRUE, datatable.print.class = TRUE, datatable.print.nrows = 50)


## constants ----

## system info

n_core <- parallel::detectCores()

## paths

dir_atlas <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"


## "metadata" / variables

## :: design info

waves <- c("wave1", "wave2")
sessions <- c("baseline", "proactive", "reactive")
sesss <- c("Bas", "Pro", "Rea")
runs <- c("run1", "run2")
#hemis <- c("L", "R")
#run_hemis <- c("run1_L", "run1_R", "run2_L", "run2_R")
#wave_dir_image <- c(wave1 = "HCP_SUBJECTS_BACKUPS", wave2 = "DMCC_Phase3")
#wave_dir_evts <- c(wave1 = "DMCC2", wave2 = "DMCC3")
colors_bias <- c("blue", "red", "purple", "white")
colors_pc50 <- c("black", "green", "pink", "yellow")
words_bias <- toupper(colors_bias)
words_pc50 <- toupper(colors_pc50)
ttypes_bias <- apply(expand.grid(colors_bias, words_bias), 1, paste0, collapse = "_")
ttypes_pc50 <- apply(expand.grid(colors_pc50, words_pc50), 1, paste0, collapse = "_")
ttypes <- c(ttypes_bias, ttypes_pc50)

n_vertex <- 20484  ## surface hcp mesh fslr5
n_tr <- c(
  baseline  = 540,
  proactive = 540,
  reactive  = 590
)  ## number of tr per subj*run for stroop task
n_trial <- c(
  baseline = 108,
  proactive = 108,
  reactive = 120
  )  ## number of trials (events) per subj*run for stroop task
n_run <- 2
n_session <- 3

## :: analytic info

glmnames <- c("lsall_1rpm", "lss_1rpm", "item_1rpm")
prewhs <- c("none", "resid", "obs", "rest")
roisets <- c("Schaefer2018_control", "Schaefer2018_network", "Schaefer2018_parcel", "Glasser2016_parcel")

core32 <- c(
  99, 127, 129, 130, 131, 132, 137, 140, 141, 142, 148, 163, 165, 182, 186, 300, 332, 333, 334, 335, 336, 337, 340, 345, 
  349, 350, 351, 352, 354, 361, 365, 387
)  ## parcel indexes in Schaefer2018




## functions ----

## misc

get_network <- function(x) gsub("^.H_(Vis|SomMot|Cont|Default|Limbic|SalVentAttn|DorsAttn)_.*", "\\1", x)
combo_paste <- function(a, b, sep = "", ...) apply(expand.grid(a, b, ...), 1, paste0, collapse = sep)
enlist <- function(nms) setNames(vector("list", length(nms)), nms)
# invert_list <- function(l) { 
#   ## https://stackoverflow.com/questions/15263146/revert-list-structure
#   ## @Josh O'Brien
#   x <- lapply(l, `[`, names(l[[1]]))  ## get sub-elements in same order
#   apply(do.call(rbind, x), 2, as.list)  ## stack and reslice
# }

Var <- function(x, dim = 1, ...) {
  ##https://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
  if(dim == 1){
     rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  } else if (dim == 2) {
     rowSums((t(x) - colMeans(x, ...))^2, ...)/(dim(x)[1] - 1)
  } else stop("Please enter valid dimension")
}


## classdef: "atlas"
## see here for more: https://github.com/mcfreund/psychomet/tree/master/in

read_atlas <- function(
    roiset = "Schaefer2018_control", 
    dir_atlas = "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"
    ) {
    
    configured <- c("Schaefer2018_control", "Schaefer2018_parcel", "Schaefer2018_network")
    if (!roiset %in% configured) stop("roiset not configured")
    
    atlas <- list()

    if (grepl("Schaefer2018", roiset)) {

        atlas$data <-
            c(
            gifti::read_gifti(file.path(dir_atlas, "Schaefer2018_400Parcels_7Networks_order_10K_L.label.gii"))$data[[1]],
            gifti::read_gifti(file.path(dir_atlas, "Schaefer2018_400Parcels_7Networks_order_10K_R.label.gii"))$data[[1]] + 200
            )
        atlas$key <- data.table::fread(here::here("in", "atlas-key_schaefer400-07.csv"))$parcel

        if (roiset == "Schaefer2018_control") {
            parcel_core32 <- atlas$key[core32]
            parcel_vis <- atlas$key[grep("_Vis_", atlas$key)]
            parcel_sommot <- atlas$key[grep("_SomMot_", atlas$key)]
            parcels_control <- c(core32 = list(parcel_core32), Vis = list(parcel_vis), SomMot = list(parcel_sommot))
            atlas$rois <- c(setNames(as.list(parcel_core32), parcel_core32), parcels_control)
        } else if ("Schaefer2018_network") {
            atlas$rois <- split(atlas$key$parcel, get_network(atlas$key$parcel))
        } else if ("Schaefer2018_parcel") {
            atlas$rois <- split(atlas$key$parcel, atlas$key$parcel)
        }

    } else if (roiset == "Glasser2016") {

    }

    atlas
    
}



## for interacting with gifti objects

extract_labels <- function(gifti) {
    stopifnot(class(gifti) == "gifti")
    vapply(gifti$data_meta, function(x) x[1, "vals"], character(1))
}

extract_data <- function(gifti, pattern = NULL) {
    stopifnot(class(gifti) == "gifti")
    data <- abind::abind(gifti$data)  ## extract
    colnames(data) <- extract_labels(gifti)
    if (!is.null(pattern)) data <- data[, grep(pattern, colnames(data))]
    data
}

concat_hemis <- function(l, hemis = c("L", "R"), pattern = NULL) {
    if (class(l) != "list" || !identical(vapply(l, class, ""), c("gifti", "gifti"))) stop("l must be list of 2 giftis")
    stopifnot(sort(hemis) == c("L", "R"))
    abind::abind(lapply(l[order(hemis)], extract_data, pattern = pattern), along = 1)
}

parcellate_data <- function(x, atlas) {
    if (!identical(class(x), "matrix")) stop("x must be of class matrix")
    if (nrow(x) != length(atlas$data)) stop("nrow(x) does not match nrow(atlas$data)")
    out <- enlist(names(atlas$rois))
    for (roi_i in seq_along(atlas$rois)) {
        which_parcels <- which(atlas$key %in% atlas$rois[[roi_i]])
        is_roi <- atlas$data %in% which_parcels
        out[[roi_i]] <- x[is_roi, ]
    }
    out
}



## filename constructors

construct_filename_gifti <- function(
    subject, wave, session, run, glmname, hemi,
    task = "Stroop", 
    base_dir = here::here("out", "glms"), 
    prefix = "STATS", 
    suffix = "REML.func.gii"
    ){
    
    arg <- as.list(environment())
    if (any(vapply(arg, length, numeric(1)) > 1)) stop("not yet configured for length>1 args")

    file.path(
        base_dir, subject, wave, "RESULTS", task, paste0(session, "_", glmname),
        paste0(prefix, "_", subject, "_", run, "_", hemi, "_", suffix)
        )
}



construct_filenames_gifti <- function(
    subjects, waves, sessions, runs, glmnames, 
    hemis = c("L", "R"), 
    task = "Stroop", 
    base_dir = here::here("out", "glms"), 
    prefix = "STATS", 
    suffix = "REML.func.gii",
    returnDT = TRUE
    ){
    
    d <- expand.grid(
        subject = subjects, wave = waves, session = sessions, run = runs, glmname = glmname,
        hemi = hemis,
        stringsAsFactors = FALSE
    )

    missing_col <- names(d)[!names(d) %in% c("subject", "wave", "session", "run", "glmname", "hemi")]
    if (length(missing_col) > 0) stop("missing column:", missing_col)
    
    filename <- vector("character", nrow(d))
    for (row_i in seq_len(nrow(d))) {
        x <- d[row_i, ]
        filename[row_i] <- construct_filename_gifti(
            subject = x$subject, wave = x$wave, 
            session = x$session, run = x$run, 
            glmname = x$glmname, hemi = x$hemi,
            task = task, base_dir = base_dir, prefix = prefix, suffix = suffix
        )
    }

    if (returnDT) {
        d$filename <- filename
        data.table::setDT(d, key = c("subject", "wave", "session", "run", "hemi"))
        d
    } else {
        filename
    }


}

construct_filename_h5 <- function(
    prefix, glmname, roiset, prewh,
    base_dir = here::here("out"), 
    suffix = ".h5"
    ){
    
    arg <- as.list(environment())
    if (any(vapply(arg, length, numeric(1)) > 1)) stop("not yet configured for length>1 args")

    file.path(base_dir, paste0(prefix, "_glm-", glmname, "_roiset-", roiset, "_prewh-", prewh, suffix))

}



## hdf5 interface functions

## writing files.
## one file per subject*wave*session*(task*run)*roi.
## dataset prefix encodes data 'type' (coefs, resids, invcov, ...).
## analysis variables (glm, prewh, ...) encoded in dataset suffix.
## a one-file-per-sample approach as discussed here: https://github.com/lgatto/MSnbase/issues/403

write_dset <- function(
    mat, dset_prefix, 
    subject, wave, session, run, roiset, roi, 
    glmname, prewh,
    base_dir = here::here("out", "parcellated", ".d")
    ) {
    id <- paste0(subject, "_", wave, "_", session, "_", run, "_", roiset, "_", roi)
    file_name <- paste0(base_dir, .Platform$file.sep, id, ".h5")
    dset_name <- paste0(dset_prefix, "_", id, "_", "_glm-", glmname, "_prewh-", prewh)
    rhdf5::h5write(mat, file = file_name, name = dset_name)
    c(file_name = file_name, dset_name = dset_name)
}

write_links <- function(names_list, base_dir = here::here("out", "parcellated")){
    master_file <- paste0(base_dir, .Platform$file.sep, "master.h5")
    if (!file.exists(master_file)) {
        fid <- rhdf5::H5Fcreate(master_file)
    } else {
        fid <- rhdf5::H5Fopen(master_file)
    }
    ## loop through files, creating a link to the dataset in each one:
    for(f_i in seq_along(names_list)) {
        f <- names_list[[f_i]]["file_name"]
        dset_name <- names_list[[f_i]]["dset_name"]
        rhdf5::H5Lcreate_external(
            target_file_name = f, 
            target_obj_name = dset_name,
            link_loc = fid,
            link_name = dset_name
            )
    }
    H5Fclose(fid)
}

