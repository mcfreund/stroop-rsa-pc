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
ttypes <- list(
    bias = sort(apply(expand.grid(colors_bias, words_bias), 1, paste0, collapse = "")),
    pc50 = sort(apply(expand.grid(colors_pc50, words_pc50), 1, paste0, collapse = ""))
)
ttypes$all <- sort(c(ttypes$bias, ttypes$pc50))
#ttypes_bias <- sort(apply(expand.grid(colors_bias, words_bias), 1, paste0, collapse = ""))
#ttypes_pc50 <- sort(apply(expand.grid(colors_pc50, words_pc50), 1, paste0, collapse = ""))
#ttypes <- sort(c(ttypes_bias, ttypes_pc50))
ttypes_by_run <- list(
    baseline = 
        list(
            run1 = sort(fread(here("in", "ttypes_baseline_run1.txt"), header = FALSE)[[1]]),
            run2 = sort(fread(here("in", "ttypes_baseline_run2.txt"), header = FALSE)[[1]])
            ),
    proactive = 
        list(
            run1 = sort(fread(here("in", "ttypes_proactive_run1.txt"), header = FALSE)[[1]]),
            run2 = sort(fread(here("in", "ttypes_proactive_run2.txt"), header = FALSE)[[1]])
            ),
    reactive = 
        list(
            run1 = sort(fread(here("in", "ttypes_reactive.txt"), header = FALSE)[[1]]),
            run2 = sort(fread(here("in", "ttypes_reactive.txt"), header = FALSE)[[1]])
            )
)
n_ttype <- c(
    baseline = 20,
    proactive = 26,
    reactive = 26
    )  ## unique ttypes PER RUN!! (gives maximum dimension of cross-run/cross-validated RDMs)
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
models <- list(
    cveuc = c("distractor", "incongruency", "target"),
    crcor = c("conjunction", "distractor", "incongruency", "target")
)
## minimum number >0 of trials per ttype (for crcor):
expected_min <- list(
    reactive_bias = 6,
    reactive_pc50 = 1,
    reactive_all = 1,
    proactive_bias = 3,
    proactive_pc50 = 3,
    proactive_all = 3,
    baseline_bias = 3,
    baseline_pc50 = 3,
    baseline_all = 3
)

## :: analytic info (expected values)

expected <- list(
    glmname = c("lsall_1rpm", "lss_1rpm", "condition_1rpm", "lsall_acompcor06_1rpm"),
    prewh   = c("none", "obsresamp", "obsresampbias", "obsresamppc50", "obsall", "obsbias", "obspc50", "obsallave"),
    roiset  = 
        c("Schaefer2018Dev", "Schaefer2018Network", "Schaefer2018Parcel", "Glasser2016Network", "Glasser2016Parcel")
)

core32 <- c(
  99, 127, 129, 130, 131, 132, 137, 140, 141, 142, 148, 163, 165, 182, 186, 300, 332, 333, 334, 335, 336, 337, 340, 345, 
  349, 350, 351, 352, 354, 361, 365, 387
)  ## parcel indexes in Schaefer2018




## functions ----

## misc

#get_network <- function(x) gsub("^.H_(Vis|SomMot|Cont|Default|Limbic|SalVentAttn|DorsAttn)_.*", "\\1", x)
# invert_list <- function(l) { 
#   ## https://stackoverflow.com/questions/15263146/revert-list-structure
#   ## @Josh O'Brien
#   x <- lapply(l, `[`, names(l[[1]]))  ## get sub-elements in same order
#   apply(do.call(rbind, x), 2, as.list)  ## stack and reslice
# }



## math/stats functions

.cvdist <- function(x1, x2, m) {
    D <- rowMeans(tcrossprod(m, x1) * tcrossprod(m, x2))  ## means to scale by num verts
    dim(D) <- sqrt(c(length(D), length(D)))  ## must be square in current implementation
    D
}

cvdist <- function(x1, x2, m = mikeutils::contrast_matrix(ncol(x1)), nms = NULL, center = FALSE, scale = FALSE) {
    if (all(dim(x1) != dim(x2))) stop("x1 and x2 must be same size")
    if (center) {
        x1 <- center(x1)
        x2 <- center(x2)
    }
    if (scale) {
        x1 <- scale2unit(x1)
        x2 <- scale2unit(x2)
    }
    D <- .cvdist(x1, x2, m)
    # attr(D, "x1_ssq") <- sqrt(colSums(x1^2))
    # attr(D, "x2_ssq") <- sqrt(colSums(x2^2))
    # attr(D, "cv_sq") <- colSums(x1 * x2)
    # attr(D, "x1_mu") <- colMeans(x1)
    # attr(D, "x2_mu") <- colMeans(x2)
    # attr(D, "n") <- nrow(x1)
    if (!is.null(nms)) dimnames(D) <- list(nms, nms)
    D
}


averaging_matrix <- function(x) {
    A <- indicator(x)
    As <- tcrossprod(A, diag(1/colSums(A)))
    colnames(As) <- unique(x)
    As
}


get_resampled_idx <- function(conditions, n_resamples, expected_min, seed = 0) {
    stopifnot(is.character(conditions) || is.numeric(n_resamples) || is.numeric(expected_min))
    set.seed(seed)

    n_conditions <- length(conditions)
    groups_list <- split(seq_along(conditions), conditions)
    resample_to <- Reduce(min, lapply(groups_list, length))
    if (resample_to != expected_min) stop("unexpected minimum n trials")
    
    replicate(n_resamples, unlist(lapply(groups_list, resample, size = resample_to)), simplify = FALSE)

}



crcor <- function(x1, x2, n_resamples, expected_min){
    stopifnot(length(dim(x1)) == 2 || length(dim(x2)) == 2)
    
    nms1 <- colnames(x1)
    nms2 <- colnames(x2)
    g1 <- unique(nms1)
    g2 <- unique(nms2)

    resample_idx1 <- 
        get_resampled_idx(conditions = nms1, n_resamples = n_resamples, expected_min = expected_min)
    resample_idx2 <- 
        get_resampled_idx(conditions = nms2, n_resamples = n_resamples, expected_min = expected_min)

    A1 <- averaging_matrix(rep(g1, each = expected_min))
    A2 <- averaging_matrix(rep(g2, each = expected_min))
    
    res <- matrix(0, ncol = ncol(A2), nrow = ncol(A1), dimnames = list(colnames(A1), colnames(A2)))
    for (ii in seq_len(n_resamples)) {
        idx1 <- resample_idx1[[ii]]
        idx2 <- resample_idx2[[ii]]
        x1i <- x1[, idx1, drop = FALSE]
        x2i <- x2[, idx2, drop = FALSE]
        res <- res + atanh(cor(x1i %*% A1, x2i %*% A2))
    }
    tanh(res / n_resamples)

}





## reading behavioral data / trial-level information

read_trialinfo <- function() {
    b <- data.table::fread(here::here("in", "behavior-and-events_stroop_2021-10-20_nice.csv"))
    dplyr::arrange(b, "subj", "wave", "session", "run", "trial_num")    ## sort to match col order of fmri beta matrix
}

## classdef: "atlas"
## see here for more: https://github.com/mcfreund/psychomet/tree/master/in

read_atlas <- function(
    roiset = "Schaefer2018Dev", 
    dir_atlas = "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"
    ) {
    
    configured <- c("Schaefer2018Dev", "Schaefer2018Parcel", "Schaefer2018Network", "Schaefer2018Parcel200")
    if (!roiset %in% configured) stop("roiset not configured")
    
    atlas <- list()

    if (grepl("Schaefer2018", roiset)) {

        if (roiset == "Schaefer2018Parcel200") {

            dir_atlas <- "/data/nil-external/ccp/freund/atlases/Schaefer200"
            L <- gifti::read_gifti(file.path(dir_atlas, "Schaefer2018_200Parcels_7Networks_order_10K_L.label.gii"))$data[[1]]
            R <- gifti::read_gifti(file.path(dir_atlas, "Schaefer2018_200Parcels_7Networks_order_10K_R.label.gii"))$data[[1]] + 100
            R[R == 100] <- 0  ## reset to zero (zero is NA...)
            atlas$data <- c(L, R)
            atlas$key <- data.table::fread(here::here("in", "atlas-key_schaefer200-07.csv"))$parcel
            atlas$rois <- split(atlas$key, atlas$key)

        } else {
            
            L <- gifti::read_gifti(file.path(dir_atlas, "Schaefer2018_400Parcels_7Networks_order_10K_L.label.gii"))$data[[1]]
            R <- gifti::read_gifti(file.path(dir_atlas, "Schaefer2018_400Parcels_7Networks_order_10K_R.label.gii"))$data[[1]] + 200
            R[R == 200] <- 0  ## reset to zero (zero is NA...)
            atlas$data <- c(L, R)
            atlas$key <- data.table::fread(here::here("in", "atlas-key_schaefer400-07.csv"))$parcel

            if (roiset == "Schaefer2018Dev") {

                parcel_core32 <- atlas$key[core32]
                parcel_vis <- atlas$key[grep("_Vis_", atlas$key)]
                parcel_sommot <- atlas$key[grep("_SomMot_", atlas$key)]
                parcels_dev <- c(core32 = list(parcel_core32), Vis = list(parcel_vis), SomMot = list(parcel_sommot))
                atlas$rois <- c(setNames(as.list(parcel_core32), parcel_core32), parcels_dev)

            } else if (roiset == "Schaefer2018Network") {

                atlas$rois <- split(atlas$key, get_network(atlas$key))

            } else if (roiset == "Schaefer2018Parcel") {

                atlas$rois <- split(atlas$key, atlas$key)

            }

        }

    } else if (roiset == "Glasser2016") {

    }

    atlas
    
}

## reading RSA models

read_model_rdm <- function(
    model, measure_type, session, ttype_subset, 
    base_dir = here::here("out", "rsa_models")
    ){
    x <- fread(
        paste0(
            base_dir, .Platform$file.sep, 
            "model_", measure_type, "_", model, "_", session, "_", ttype_subset, ".csv")
        )
    R <- as.matrix(x[, -1])
    rownames(R) <- x$run1
    R
}

read_model_xmat <- function(
    measure, session, ttype_subset, 
    base_dir = here::here("out", "rsa_models"), 
    .models = models[[measure]]
    ) {
    mods <- enlist(.models)
    for (mod in .models) {
        mods[[mod]] <- read_model_rdm(
            model = mod, 
            measure_type = switch(measure, crcor = "similarity", cveuc = "cvdistance"),
            session = session,
            ttype_subset = ttype_subset
            )
    }
    cbind(intercept = 1, do.call(cbind, lapply(mods, switch(measure, crcor = c, cveuc = mfutils::squareform))))
}

## tidying regression output

tidy_model <- function(B, terms, outcomes) {
    B <- as.data.table(B)
    names(B) <- outcomes
    B$term <- terms
    B
 }


## filename constructors

construct_filename_gifti <- function(
    subject, wave, session, run, glmname, hemi,
    task = "Stroop", 
    base_dir = here::here("out", "glms"), 
    prefix = "STATS", 
    suffix = "_REML.func.gii"
    ){
    
    arg <- as.list(environment())
    if (any(vapply(arg, length, numeric(1)) > 1)) stop("not yet configured for length>1 args")

    file.path(
        base_dir, subject, wave, "RESULTS", task, paste0(session, "_", glmname),
        paste0(prefix, "_", subject, "_", run, "_", hemi, suffix)
        )
}



construct_filenames_gifti <- function(
    subjects, waves, sessions, runs, glmnames, 
    hemis = c("L", "R"), 
    task = "Stroop", 
    base_dir = here::here("out", "glms"), 
    prefix = "STATS", 
    suffix = "_REML.func.gii",
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


construct_filename_weights <- function(
    measure, subjlist, glmname, ttype_subset, roiset, prewh,
    prefix = "weights", base_dir = here::here("out", "res"), suffix = ""
    ) {
    paste0(
        base_dir, .Platform$file.sep, 
        prefix, "-", measure, "__subjlist-", subjlist, "__glm-", glmname, "__ttype-", ttype_subset,
        "__roiset-", roiset, "__prewh-", prewh, suffix, ".csv"
        )
}



## hdf5 interface functions

## reading/writing parcellated data files (matrices built from segmented gifti images).

## dataset_prefix encodes data 'type' (coefs, resids, invcov, ...).
## analysis variables (glm, prewh, ...) encoded in dataset_suffix.
## these functions are structured to support a one file per "sample" approach: subject*wave*session*(task*run)*roi
## due to separate files, the separate "samples" can be read/written in parallel.
## (see here for helpful discussion of issues: https://github.com/lgatto/MSnbase/issues/403)
## in the present case, this allows essentially all central operations of a script to be wrapped in a parallelizing 
## function (such as foreach::foreach()), simplifying development of fast code.
## the downside of this approach is that the directory containing these files will be very wide, with many thousands of
## files.
## thus this directory (parcellated/.d) will be hidden.
## to aid in accessing (subsets of) these files/datasets, a master.h5 file containing external links to each dataset
## will be written in the top directory (parcellated).
## importantly, this master file will need to be created *outside* of the parallelizing loop (in serial).
## for convenience, write_dset() returns the filename and dataset name of the written data as a named character vector.
## write_links() therefore a list of such character vectors, loops over them serially, and creates this master file.
## read_master() reads this master file of links and does some string editing to make subsetting easier.
## an additional option includes saving of colnames.
## write_dset() optionally saves colnames of matrix as dataset attribute.
## read_dset() optionally reads colnames and adds as dimnames attribute to R matrix.
## read_dset() uses external links within the master file to access the file/dataset.

# .construct_filename_rdm
# .construct_filename_gifti
# .construct_filename_coefs


# .construct_dsetname_coefs
# .construct_dsetname_rdm

construct_filename_rdm <- function(
    measure, glmname, ttype_subset, roiset, prewh, 
    base_dir = here::here("out", "res")
    ){
    paste0(
        base_dir, .Platform$file.sep, 
        "rdm-", measure, "__glm-", glmname, "__ttype-", ttype_subset,  "__roiset-", roiset, "__prewh-", prewh, ".h5"
        )
}



construct_filename_h5 <- function(
    subject, wave, session, run, roiset, roi, 
    base_dir = here::here("out", "parcellated", ".d")
    ){
    paste0(base_dir, .Platform$file.sep, subject, "_", wave, "_", session, "_", run, "_", roiset, "_", roi, ".h5")
}
construct_dsetname_h5 <- function(
    dset_prefix, subject, wave, session, run, roiset, roi, 
    glmname, prewh,
    base_dir = here::here("out", "parcellated", ".d")
    ){
    paste0(
        dset_prefix, "__", subject, "__", wave, "__", session, "__", run, "__", roiset, "__", roi, "__", 
        "glm-", glmname, "__prewh-", prewh
        )
}


.write_dset_dimnames <- function(mat, file_name, dset_name, idx) {
    fid <- rhdf5::H5Fopen(file_name)
    did <- rhdf5::H5Dopen(fid, dset_name)
    rhdf5::h5writeAttribute(dimnames(mat)[[idx]], did, switch(idx, "rownames", "colnames"))
    rhdf5::H5Dclose(did)
    rhdf5::H5Fclose(fid)
}

write_dset <- function(
    mat, dset_prefix, subject, wave, session, run, roiset, roi, glmname, prewh,
    base_dir = here::here("out", "parcellated", ".d"),
    write_colnames = FALSE,
    write_rownames = FALSE
    ) {
    file_name <- construct_filename_h5(subject, wave, session, run, roiset, roi, base_dir)
    dset_name <- construct_dsetname_h5(dset_prefix, subject, wave, session, run, roiset, roi, glmname, prewh, base_dir)
    rhdf5::h5write(mat, file = file_name, name = dset_name)
    if (write_rownames) .write_dset_dimnames(mat, file_name, dset_name, 1)
    if (write_colnames) .write_dset_dimnames(mat, file_name, dset_name, 2)
    c(file_name = file_name, dset_name = dset_name)
}

construct_filenames_h5 <- function(
    prefix, subjects, waves, sessions, rois, runs, glmname, prewh
    ){    
    input <- expand.grid(
        subj = subjects, wave = waves, session = sessions, roi = rois, run = runs, 
        stringsAsFactors = FALSE
        )
    setDT(input)
    input[, file_name := construct_filename_h5(subj, wave, session, run, roiset, roi)]
    input[, dset_name := construct_dsetname_h5(prefix, subj, wave, session, run, roiset, roi, glmname, prewh)]
    input
}


construct_filename_datainfo <- function(
    prefix, subjlist, glmname, roiset, prewh, base_dir = here::here("out", "parcellated")
    ) {
    paste0(
        base_dir, .Platform$file.sep, "datainfo-", prefix, "__subjlist-", subjlist, "__glm-", glmname,
        "__roiset-", roiset, "__prewh-", prewh, ".csv"
    )
}


parse_dset_name <- function(
    x, nms = c("prefix", "subj", "wave", "session", "run", "roiset", "roi", "glm", "prewh")
    ) {
    tidyr::separate(x, name, nms, sep = "__")
}

read_dset <- function(file_name, dset_name, read_colnames = TRUE) {
    mat <- rhdf5::h5read(file_name, dset_name, read.attributes = TRUE)
    if (read_colnames) {
        colnames(mat) <- attr(mat, "colnames")
        attr(mat, "colnames") <- NULL
    }
    mat
}


read_rdms <- function(
    .measure, .glmname, .roiset, .prewh, .subjects, .session, .waves,
    .ttype_subset = ttype_subset,
    .ttypes1 = NULL, 
    .ttypes2 = NULL, 
    .rois = names(atlas$rois)
    ) {
    stopifnot(length(.ttypes1) == length(.ttypes2))
    on.exit(rhdf5::h5closeAll())

    n_ttype <- length(.ttypes1)
    n_sub <- length(.subjects)
    n_wav <- length(.waves)
    n_roi <- length(.rois)

    A <- array(
        NA, 
        dim = c(n_ttype, n_ttype, n_roi, n_sub, n_wav),
        dimnames = list(dim1 = .ttypes1, dim2 = .ttypes2, roi = .rois, subject = .subjects, wave = .waves)
        )
    
    fname <- construct_filename_rdm(
        measure = .measure, glmname = .glmname, ttype_subset = .ttype_subset, roiset = .roiset, prewh = .prewh
        )
    fid <- rhdf5::H5Fopen(fname)
    for (.subject in .subjects) {
        for (.wave in .waves) {
            dset_name <- paste0(.subject, "__", .wave, "__", .session)
            did <- rhdf5::H5Dopen(fid, dset_name)
            A[, , , .subject, .wave] <- rhdf5::H5Dread(did)
            rhdf5::H5Dclose(did)
        }
    }
    

    A

}



## resampling functions



resample_apply_combine <- function(
    x, resample_idx,
     apply_fun = identity, 
     combine_fun = function(.x) Reduce("+", .x) / length(.x),
     outdim = NULL
     ) {
    stopifnot(is.matrix(x) || is.matrix(resample_idx))

    n_resamples <- nrow(resample_idx)

    if (is.function(combine_fun)) {
        res <- vector("list", n_resamples)
        for (ii in seq_len(n_resamples)) {
            idx <- resample_idx[ii, ]
            xii <- x[, idx, drop = FALSE]
            res[[ii]] <- apply_fun(xii)
        }
        res <- combine_fun(res)
    } else if (combine_fun == "iterative_add") {
        if (is.null(outdim)) stop("outdim must be specified when iterative_add == TRUE")
        res <- list(matrix(0, nrow = outdim[1], ncol = outdim[2]))  ## wrap in list so modify in place
        #res <- matrix(0, nrow = outdim[1], ncol = outdim[2])
        for (ii in seq_len(n_resamples)) {
            idx <- resample_idx[ii, ]
            xii <- x[, idx, drop = FALSE]
            res[[1]] <- res[[1]] + apply_fun(xii)/n_resamples
            #res <- res + apply_fun(xii)/n_resamples
        }
        res <- res[[1]]
    }

    res

}

