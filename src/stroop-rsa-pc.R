## 
## defines functions for use in this project

## class def
## https://adv-r.hadley.nz/s3.html, sec 13.3.1

new_parcellated_list  <- function(x = list(), glm_name = character(), roi_set = character(), prewhitened = logical(), shrinkage = numeric()) {
    stopifnot(is.list(x) || is.character(glm_name) || is.character(roi_set) || is.logical(prewhitened) || is.numeric(shrinkage))
    structure(x, class = "parcellated_list", glm_name = glm_name, roi_set = roi_set, prewhitened = prewhitened, shrinkage = shrinkage)
}


x <- parc_gii
validate_parcellated_list <- function(x) {
    ## check types within list
    ## check consistency of features across folds, obs across rois
    ## check values of attr
    
    if (!identical("list", unique(vapply(x, class, character(1))))) stop("x is not nested list")
    if (any(!c("data", "labels") %in% names(x))) stop("x does not contain 'data' or 'labels' names")
    
    x_data <- x$data
    x_labels <- x$labels
    if (!identical("list", unique(vapply(x_data, class, character(1))))) stop("x$data is not nested list")
    if (!identical("character", unique(vapply(x_labels, class, character(1))))) stop("x$labels is not list of character vectors")
    
    n_fold_data <- unique(vapply(x_data, length, numeric(1)))
    n_fold_label <- length(x_labels)
    if (length(n_fold_data) != 1L) stop("Contains differing numbers of folds across ROIs.")
    if (n_fold_data != n_fold_label) stop("Mismached folds in data and labels.")

    if (any(duplicated(names(x_data)))) stop("x$data contains duplicate ROI names")

    u <- unlist(x_data, recursive = FALSE)
    if (any(duplicated(u))) stop("x$data contains duplicate roi*folds")
    vapply(u, class, character(1))
    n_obs <- unique(vapply(u, ncol, numeric(1)))
    if (length(n_obs) != 1L) stop("x$data number of observations")
    
    





}

invert_list <- function(l) { 
  ## https://stackoverflow.com/questions/15263146/revert-list-structure
  ## @Josh O'Brien
  x <- lapply(l, `[`, names(l[[1]]))  ## get sub-elements in same order
  apply(do.call(rbind, x), 2, as.list)  ## stack and reslice
}

as_parcellated_list  <- function(giftis, fold_hemis, atlas, labels_from_data = TRUE, pattern = NULL) {
    ## giftis: list of gifti images of length N_folds*N_hemi (one element per fold*hemi)
    ## hemis, folds: vectors of length length(giftis), specifying hemi and fold value per element of gifti
    ## atlas: list of data and key of parcellation atlas
   
   fold_hemis <- strsplit(fold_hemis, "_")
   folds <- unlist(lapply(fold_hemis, "[", 1))
   hemis <- unlist(lapply(fold_hemis, "[", 2))

    if (!(length(giftis) == length(folds) && length(folds) == length(hemis))) stop("giftis must have length equal to length(folds) == length(hemis)")
    if (!(class(giftis[[1]]) == class(giftis[[2]]) && class(giftis[[1]]) == "gifti")) stop("elements of giftis must be of class gifti")
    #if (isFALSE(all.equal(giftis[[1]]$data_info$Dim0, giftis[[2]]$data_info$Dim0))) stop("giftis must have same dimensions")
    if (!is.character(folds)) stop("folds must be class character")
    if (!is.character(hemis)) stop("hemis must be class character")
    if (!any(hemis %in% c("L", "R"))) stop("hemis not in L, R")
    if (length(hemis) %% 2 != 0L) stop("length(hemis) must be even number")
    if (!is.logical(labels_from_data)) stop("labels_from_data must be logical")
    if (!(is.character(pattern) || is.null(pattern))) stop("pattern must be character or NULL")

    ## extract data, labels

    l <- split(giftis, folds)
    hemis_split <- split(hemis, folds)  ## for reordering 
    unique_folds <- unique(folds)
    d <- enlist(unique_folds)
    labels <- enlist(unique_folds)
    for (fold_i in seq_along(l)) {

        l_val <- lapply(l[[fold_i]], function(x) abind::abind(x$data))  ## extract
        l_val <- l_val[order(hemis_split[[fold_i]])]  ## ensure L then R
        l_val <- abind::abind(l_val, along = 1)  ## L then R concat
        
        if (labels_from_data) {
            labs <- vapply(l[[fold_i]][[1]]$data_meta, function(x) x[1, "vals"], character(1))
            #labs_R <- vapply(l[[fold_i]][[2]]$data_meta, function(x) x[1, "vals"], character(1))  ## add to is.parcellated_image()
            #if (!identical(labs, labs_R)) stop("labs not identical across hemisphere")
            colnames(l_val) <- labs
            if (!is.null(pattern)) {
                idx <- grep(pattern, labs)
                l_val <- l_val[, idx]
                labs <- labs[idx]
            }
            labels[[fold_i]] <- labs
        }

        d[[fold_i]] <- l_val

    }
    
    ## parcellate data
    
    data <- enlist(names(atlas$rois))
    for (roi_i in seq_along(atlas$rois)) {
        which_parcels <- which(atlas$key %in% atlas$rois[[roi_i]])
        is_roi <- atlas$data %in% which_parcels
        data[[roi_i]] <- lapply(d, function(x) x[is_roi, ])
    }
    
    ## get "good list"
    
    get_good_vertices <- function(x) {
        n_vertex_roi <- unique(vapply(x, nrow, numeric(1)))
        if (length(n_vertex_roi) != 1) stop ("something very wrong")
        is_good_vertex_folds <- vapply(x, function(y) apply(y, 1, var) > 0, logical(n_vertex_roi))
        rowSums(is_good_vertex_folds) > 0  ## collapse across folds (must be good in all)
    }
    good_vertices <- lapply(data, get_good_vertices)

    out <- list(data = data, good_vertices = good_vertices, labels = labels)
    class(out) <- c("parcellated_list", "list")

    out

}

is.parcellated_list <- function(x) inherits(x, "parcellated_list")

print.parcellated_list <- function(x) {

    n_roi <- length(x$data)
    n_fold <- length(x$data[[1]])
    dims <- vapply(x$data, function(x) dim(x[[1]]), numeric(2))
    range_vertex <- range(dims[1, ])
    n_trialtype <- unique(dims[2, ])
    # cat(length(x$data[[1]]), "folds \n")
    # cat(length(x$data[[1]]), "folds \n")
    out_string <- paste0(
        "parcellated list of ", n_roi, " ROIs (", min(range_vertex), "-", max(range_vertex), " features by ", n_trialtype, " conditions), ", 
        "each with ", n_fold, " folds\n"
    )

    cat(out_string)

}

coef.parcellated_list <- function(x, rois = NULL, folds = NULL, good_vertices = TRUE, reduce = TRUE) {
    
    get_good_vertices_loop <- function(.data_roi, .good_vertices_roi) lapply(.data_roi, function(.x) .x[.good_vertices_roi, ])

    if (is.null(rois)) {
        out <- x$data
    } else {
        out <- x$data[rois]
    }

    if (!is.null(folds)) out <- lapply(out, function(a) a[folds])

    if (good_vertices) {
        if (!is.null(rois)) {
            good_vertices_rois <- x$good_vertices[rois, drop = FALSE]
        } else {
            good_vertices_rois <- x$good_vertices
        }
        out <- Map(function(a, b) get_good_vertices_loop(a, b), out, good_vertices_rois)
    }

    if (reduce) {
        if (length(rois) == 1) {  ## when only one ROI selected
            out <- out[[1]]
            if (length(folds) == 1) out <- out[[1]]
        } else if (length(folds) == 1) {
            out <- lapply(out, function(a) a[[1]])
        }
    }

    return(out)

}



relabel <- function(x, nms, rename_cols = FALSE) {
    if (!is.list(nms)) stop("nms must be list")
    if (length(x$labels) != length(nms)) stop("n_folds in new labels must equal n_folds in old labels")
    if (!identical(lapply(x$labels, length), lapply(nms, length))) stop("n_folds and names in new labels must equal those in old labels")  ## ensure check to ncol in is.parcellated_image()


    x$labels <- nms

    if (rename_cols) {
        for (roi_i in seq_along(x$data)) {
            for (fold_i in seq_along(x$data[[roi_i]])) colnames(x$data[[roi_i]][[fold_i]]) <- nms[[fold_i]]
        }
    }

    x

}


parcellated_lapply <- function(x, f, ..., as_parcellated_list = TRUE) {
    if (!is.parcellated_list(x)) stop("not parcellated_list")  ## add check to is....() for same number of cols across ROIs.

    u <- unlist(x$data, recursive = FALSE)
    res <- parallel::mclapply(u, f, ...)
    names(res) <- unlist(lapply(x$data, names), use.names = FALSE)
    x$data <- split(res, names(x$data))

    if (as_parcellated_list && !is.parcellated_list(x)) stop("result no longer parcellated_list")
    
    x

}

as.array.parcellated_list <- function(x) {
    if (!is.parcellated_list(x)) stop("not parcellated_list")  ## add check to is....() for same number of cols across ROIs.
    #n_cols <- unlist(lapply(x$data, function(x) lapply(x, ncol)))
    #n_unique_ncols <- length(unique(n_cols))
    n_rows <- unlist(lapply(x$data, function(x) lapply(x, nrow)))
    n_unique_nrows <- length(unique(n_rows))
    if (n_unique_nrows != 1L) stop("n_rows must be matched across folds")

    nms <- list(
        Var1 = rownames(x$data[[1]][[1]]),
        Var2 = colnames(x$data[[1]][[1]]),
        fold = names(x$data[[1]]),
        roi = names(x$data)
    )
    
    u <- unlist(x$data, recursive = FALSE)
    u <- u[sort(combo_paste(nms$roi, nms$fold, sep = "."))]  ## sort b/c unlist screws up order?

    a <- abind::abind(u, along = 0)
    ap <- aperm(a, c(2, 3, 1))
    dim(ap) <- lapply(nms, length)
    dimnames(ap) <- nms
    
    ap

}


## read atlas

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

## misc




get_network <- function(x) gsub("^.H_(Vis|SomMot|Cont|Default|Limbic|SalVentAttn|DorsAttn)_.*", "\\1", x)
combo_paste <- function(a, b, sep = "", ...) apply(expand.grid(a, b, ...), 1, paste0, collapse = sep)
enlist <- function(nms) setNames(vector("list", length(nms)), nms)
