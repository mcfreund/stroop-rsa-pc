## 
## defines functions for use in this project

## class def


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
    
    data <- enlist(names(rois))
    for (roi_i in seq_along(rois)) {
        which_parcels <- which(atlas$key$parcel %in% rois[[roi_i]])
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
