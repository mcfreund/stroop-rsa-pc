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



## classef: "atlas"
## see here for more: https://github.com/mcfreund/psychomet/tree/master/in

read_atlas <- function(roiset = "Schaefer2018_control", dir_atlas = "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/") {
    
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
    task = "Stroop", base_dir = here::here("out", "glms"), prefix = "STATS", suffix = "REML.func.gii"
    ){
    
    arg <- as.list(environment())
    if (any(vapply(arg, length, numeric(1)) > 1)) stop("not yet configured for length>1 args")

    file.path(
        base_dir, subject, wave, "RESULTS", task, paste0(session, "_", glmname),
        paste0(prefix, "_", subject, "_", run, "_", hemi, "_", suffix)
        )
}



construct_filenames_gifti <- function(
    d,
    # subjects, waves, sessions, runs, glmnames, 
    # hemis = c("L", "R"), 
    task = "Stroop", base_dir = here::here("out", "glms"), prefix = "STATS", suffix = "REML.func.gii"
    ){
    
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

    filename

}

construct_filename_h5 <- function(
    prefix, glmname, roiset, prewh,
    base_dir = here::here("out"), suffix = ".h5"
    ){
    
    arg <- as.list(environment())
    if (any(vapply(arg, length, numeric(1)) > 1)) stop("not yet configured for length>1 args")

    file.path(base_dir, paste0(prefix, "_glm-", glmname, "_roiset-", roiset, "_prewh-", prewh, suffix))

}



## I/O with hdf5 (wrappers for rhdf5 functions)

create_h5groups <- function(filename, subjects, waves = NULL, sessions = NULL, runs = NULL, rois = NULL){
    
    f <- function(filename, ...) {
        x <- expand.grid(..., stringsAsFactors = FALSE)
        paths <- paste0("/", apply(x, 1, paste0, collapse = "/"))
        lapply(paths, rhdf5::h5createGroup, file = filename)
    }

    f(filename, subjects)
    if (!is.null(waves)) {
      f(filename, subjects, waves)
        if (!is.null(sessions)) {
            f(filename, subjects, waves, sessions)
            if (!is.null(runs)) {
                f(filename, subjects, waves, sessions, runs)
                if (!is.null(rois)) {
                    f(filename, subjects, waves, sessions, runs, rois)
                }
            }
        }
    }
    
    rhdf5::h5closeAll()

} 







# write_nested_h5 <- function(x, colname_data, file_name, data_name = colname_data) {
    
#     d <- x$data
#     files <- apply(d[, .(subject, session, wave, roi, fold)], 1, paste0, collapse = "/")
#     data <- d[[colname_data]]
    
#     UseMethod("write_nested_h5")
# }

# write_nested_h5.parcellated_data <- function(x, colname_data, file_name, data_name = colname_data) {

#     for (file_i in seq_along(files)) {
#         g <- files[file_i]
#         x <- data[[file_i]]
#         rhdf5::h5write(x, file_name, paste0(c(g, data_name), collapse = "/"))
#     }

# }


# write_nested_h5.parcellated_image <- function(x, colname_data, file_name, data_name = colname_data) {

#     bad_vertices <- d$bad_vertices

#     for (file_i in seq_along(files)) {
#         g <- files[file_i]
#         x <- data[[file_i]]
#         if (length(bad_vertices[[file_i]]) > 0) {
#             rhdf5::h5write(bad_vertices[[file_i]], file_name, paste0(c(g, "bad_vertices"), collapse = "/"))
#         }
#         rhdf5::h5write(x, file_name, paste0(c(g, data_name), collapse = "/"))
#     }

# }



# construct_filename <- function(
#     prefix,
#     filetype,
#     pdata = NULL,
#     data_type = NULL,
#     glm_name = NULL,
#     roi_set = NULL,
#     prewhitened = NULL,
#     only_good_verts = NULL
#     ) {

#     strings <- c(data_type, glm_name, roi_set, prewhitened)
    
#     if (!is.null(pdata)) {

#         if (any(!is.null(strings))) stop("strings and pdata must not be both specified")
#         glm_name <- pdata$metadata$glm_name
#         roi_set <- pdata$metadata$roi_set
#         prewhitened <- pdata$metadata$prewhitened
#         shrinkage_var <- pdata$metadata$shrinkage_var
#         shrinkage_cov <- pdata$metadata$shrinkage_cov
#         if (is.parcellated_image(pdata)) {
#             prefix <- paste0("pimage-", prefix)
#             bad_vertices <- ""
#         } else if (is.parcellated_data(pdata)) {
#             bad_vertices <- paste0("badvertices-", switch(only_good_verts + 1, "incl", "rm"))
#             prefix <- paste0("pdata-", prefix)
#         } else stop("if pdata supplied must be either parcellated_image or parcellated_data")

#     } else {
        
#         if (any(is.null(strings))) stop("missing character arg")
        
#         expected <- list(
#             data_type = c("pdata", "pimage"),
#             glm_name = c("lsall_1rpm", "lss_1rpm", "item_1rpm"),
#             roi_set = c("Schaefer2018_control", "Schaefer2018_network", "Schaefer2018_parcel"),
#             prewhitened = c("none", "resid", "obs"),
#             only_good_verts = c(TRUE, FALSE)
#         )

#         if (!data_type %in% expected$data_type) stop("data_type must be one of ", paste0(expected$data_type, sep = " "))
#         if (!glm_name %in% expected$glm_name) stop("data_type must be one of ", paste0(expected$glm_name, sep = " "))
#         if (!roi_set %in% expected$roi_set) stop("data_type must be one of ", paste0(expected$roi_set, sep = " "))
#         if (!prewhitened %in% expected$prewhitened) stop("data_type must be one of ", paste0(expected$prewhitened, sep = " "))
#         if (data_type == "pdata" && is.null(only_good_verts)) stop("must specify only_good_verts with data_type = pdata")
#         prefix <- paste0(data_type, "-", prefix)
#         if (data_type == "pdata") {
#             bad_vertices <- paste0("badvertices-", switch(only_good_verts + 1, "incl", "rm"))
#         } else bad_vertices <- ""
#     }

#     paste0(
#         prefix, "_", 
#         "glm-", glm_name, 
#         "_roiset-", roi_set, 
#         "_prewh-", prewhitened,
#         bad_vertices,
#         filetype
#         )

# }




## classdef: "parcellated_data"
## https://adv-r.hadley.nz/s3.html, sec 13.3.1

## parellated_data
## list of lenth two: data and labels
## data: a data.table 
##      - contains a list-column of feature (e.g., vertex) by obs (e.g., trial) matrices, one matrix per row.
##      - other columns specify information about each row's data matrix. e.g., roi, fold, subject, session, wave, ...
##      - in all likelihood, other columns contain infomation that could vary across rows.
## metadata: a list of infor that likely does not vary across rows of data
##      - e.g., glm_name, prewhitened, shrinkage_var, shrinkage_cov, ...
##
## parcellated_image
## parcellated_data object whos features are vertices and that contains an additional list-column of "bad_verts".
## bad_verts are integer vectors that index vertices with no measured fMRI signal.

# new_parcellated_data <- function(
#     x = data.table::data.table(),
#     glm_name = character(), 
#     roi_set = character(), 
#     prewhitened = character(), 
#     shrinkage_var = numeric(),
#     shrinkage_cov = numeric()
#     ) {
#     stopifnot(
#         data.table::is.data.table(x) && is.character(glm_name) && is.character(roi_set) && is.character(prewhitened) && 
#         is.numeric(shrinkage_var) && is.numeric(shrinkage_cov)
#         )
#     obj <- list(
#         data = x,
#         metadata = list(glm_name = glm_name, roi_set = roi_set, prewhitened = prewhitened, shrinkage_var = shrinkage_var, shrinkage_cov = shrinkage_cov)
#         )
#     structure(obj, class = c("parcellated_data", "list"))
# }

# new_parcellated_image <- function(...) {
#     x <- new_parcellated_data(...)
#     class(x) <- c("parcellated_image", class(x))
#     x
# }

#validate_parcellated_data <- function()  ## TODO: ensure all attributes are of length 1
#validate_parcellated_image  ## check bad_vertices


# validate_parcellated_image <- function(x) {
#     ## check values within data.table
#     ## check consistency of features across folds, obs across rois
#     ## check values of attr
    
#     if (!identical("list", unique(vapply(x, class, character(1))))) stop("x is not nested list")
#     if (any(!c("data", "labels") %in% names(x))) stop("x does not contain 'data' or 'labels' names")
    
#     x_data <- x$data
#     x_labels <- x$labels
#     if (!identical("list", unique(vapply(x_data, class, character(1))))) stop("x$data is not nested list")
#     if (!identical("character", unique(vapply(x_labels, class, character(1))))) stop("x$labels is not list of character vectors")
    
#     n_fold_data <- unique(vapply(x_data, length, numeric(1)))
#     n_fold_label <- length(x_labels)
#     if (length(n_fold_data) != 1L) stop("Contains differing numbers of folds across ROIs.")
#     if (n_fold_data != n_fold_label) stop("Mismached folds in data and labels.")

#     if (any(duplicated(names(x_data)))) stop("x$data contains duplicate ROI names")

#     u <- unlist(x_data, recursive = FALSE)
#     if (any(duplicated(u))) stop("x$data contains duplicate roi*folds")
#     vapply(u, class, character(1))
#     n_obs <- unique(vapply(u, ncol, numeric(1)))
#     if (length(n_obs) != 1L) stop("x$data number of observations")
    
# }

#merge.parcellated_data


## converting list of giftis to parcellated_image object

giftis_to_parcellated_image <- function(
    giftis, folds, hemis, atlas, colname_data, 
    subject, wave, session, task,
    glm_name, roi_set, prewhitened, shrinkage_var, shrinkage_cov,
    pattern = NULL
    ) {
    ## giftis: list of gifti images of length N_folds*N_hemi (one element per fold*hemi)
    ## hemis, folds: vectors of length length(giftis), specifying hemi and fold value per element of gifti
    ## atlas: list of data and key of parcellation atlas
   
    ## type checks
   
    if (!(length(giftis) == length(folds) && length(folds) == length(hemis))) stop("giftis must have length equal to length(folds) == length(hemis)")
    if (class(giftis) != "list") stop("giftis must be of class list")
    if (!identical(unique(vapply(giftis, class, "")), "gifti")) stop("giftis must contain only elements of class gifti")
    if (!is.character(folds)) stop("folds must be class character")
    if (!is.character(hemis)) stop("hemis must be class character")
    if (!any(hemis %in% c("L", "R"))) stop("hemis not in L, R")
    if (length(hemis) %% 2 != 0L) stop("length(hemis) must be even number")
    ## TODO: ADD ATLAS CHECKS
    if (!(is.character(pattern) || is.null(pattern))) stop("pattern must be character or NULL")
    
    ## extract data and labels, concatenate hemispheres, and parcellate
    l <- split(giftis, folds)
    hemis_split <- split(hemis, folds)  ## for reordering
    data <- Map(function(a, b, pattern) concat_hemis(a, b, pattern), a = l, b = hemis_split, pattern = pattern)
    parcs <- lapply(data, parcellate_data, atlas = atlas)
    
    ## combine in data.table

    u <- unlist(parcs, recursive = FALSE)
    dt <- data.table::data.table(  ## get roi/fold strings from names(u) to ensure order match w/ u
        roi = gsub(paste0(unique(folds), "\\.", collapse = "|"), "", names(u)),
        fold = gsub(paste0("\\.", unique(names(atlas$rois)), collapse = "|"), "", names(u)),
        subject = subject,
        wave = wave,
        session = session,
        task = task
    )
    data.table::setkey(dt, roi, fold, subject, wave, session, task)
    #names(u)[match(paste0(dt$fold, ".", dt$roi), names(u))] == paste0(dt$fold, ".", dt$roi)
    u <- u[match(paste0(dt$fold, ".", dt$roi), names(u))]  ## reorder u to match dt
    data.table::set(dt, j = colname_data, value = u)  ## add u
    
    ## "bad list"

    bad_vertices <- lapply(u, function(.x) which(!apply(.x, 1, var) > 0))  ## verts with no variance over obs in each fold
    set(dt, j = "bad_vertices", value = bad_vertices)
    #bad_vertices <- lapply(split(bad_vertices, dt$roi), function(x) Reduce(union, x))  ## union over folds
    #bad_vertices <- bad_vertices[lapply(bad_vertices, length) > 0]

    new_parcellated_image(
        dt,
        glm_name = glm_name, roi_set = roi_set, prewhitened = prewhitened, shrinkage_var = shrinkage_var, shrinkage_cov = shrinkage_cov
        )

}

is.parcellated_data <- function(x) inherits(x, "parcellated_data")
is.parcellated_image <- function(x) inherits(x, "parcellated_image")



## methods

print.parcellated_data <- function(x) print(x$data)
print.parcellated_image <- function(x) {
    print(x$data)
    cat( 
        sum(vapply(x$data$bad_vertices, length, numeric(1))), "bad vertices in",  
        length(x$data$bad_vertices), "ROI(s)\n"
        )
}

## coef.parcellated_image
## coef.parcellated_data
## as_array
## lapply_parcellated_data
## Map_parcellated_data
## relabel

# parcellated_lapply <- function(x, f, ..., as_parcellated_list = TRUE) {
#     if (!is.parcellated_data(x)) stop("not parcellated_data")  ## add check to is....() for same number of cols across ROIs.
#     x$data <- parallel::mclapply(x$data, f, ...)
#     if (as_parcellated_list && !is.parcellated_list(x)) stop("result no longer parcellated_list")    
#     x
# }

# coef.parcellated_list <- function(x, rois = NULL, folds = NULL, good_vertices = TRUE, reduce = TRUE) {
    
#     get_good_vertices_loop <- function(.data_roi, .good_vertices_roi) lapply(.data_roi, function(.x) .x[.good_vertices_roi, ])

#     if (is.null(rois)) {
#         out <- x$data
#     } else {
#         out <- x$data[rois]
#     }

#     if (!is.null(folds)) out <- lapply(out, function(a) a[folds])

#     if (good_vertices) {
#         if (!is.null(rois)) {
#             good_vertices_rois <- x$good_vertices[rois, drop = FALSE]
#         } else {
#             good_vertices_rois <- x$good_vertices
#         }
#         out <- Map(function(a, b) get_good_vertices_loop(a, b), out, good_vertices_rois)
#     }

#     if (reduce) {
#         if (length(rois) == 1) {  ## when only one ROI selected
#             out <- out[[1]]
#             if (length(folds) == 1) out <- out[[1]]
#         } else if (length(folds) == 1) {
#             out <- lapply(out, function(a) a[[1]])
#         }
#     }

#     return(out)

# }


relabel <- function(pdata, colname_data, nms) {
    if (!is.character(colname_data) && length(colname_data) == 1L) stop("colname must be character(1)")
    if (!is.list(nms)) stop("nms must be list")
    if (length(nms) != length(unique(pdata$data$fold))) stop("length(nms) must equal number of folds")
    
    for (row_i in seq_len(nrow(pdata$data))) {
        if (ncol(pdata$data[[colname_data]][row_i][[1]]) != length(nms[[pdata$data[row_i, fold]]])) stop(paste0("ncol mismatch: row ", row_i))
        colnames(pdata$data[[colname_data]][row_i][[1]]) <- nms[[pdata$data[row_i, fold]]]
    }
    
    pdata

}





# as.array.parcellated_list <- function(x) {
#     if (!is.parcellated_list(x)) stop("not parcellated_list")  ## add check to is....() for same number of cols across ROIs.
#     #n_cols <- unlist(lapply(x$data, function(x) lapply(x, ncol)))
#     #n_unique_ncols <- length(unique(n_cols))
#     n_rows <- unlist(lapply(x$data, function(x) lapply(x, nrow)))
#     n_unique_nrows <- length(unique(n_rows))
#     if (n_unique_nrows != 1L) stop("n_rows must be matched across folds")

#     nms <- list(
#         Var1 = rownames(x$data[[1]][[1]]),
#         Var2 = colnames(x$data[[1]][[1]]),
#         fold = names(x$data[[1]]),
#         roi = names(x$data)
#     )
    
#     u <- unlist(x$data, recursive = FALSE)
#     u <- u[sort(combo_paste(nms$roi, nms$fold, sep = "."))]  ## sort b/c unlist screws up order?

#     a <- abind::abind(u, along = 0)
#     ap <- aperm(a, c(2, 3, 1))
#     dim(ap) <- lapply(nms, length)
#     dimnames(ap) <- nms
    
#     ap

# }

