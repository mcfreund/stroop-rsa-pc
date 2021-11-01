## 
## defines functions for use in this project

## class def

# pattern <- "_Coef"
# labels_from_data <- TRUE
# giftis <- gii

parcellated_image  <- function(giftis, folds, hemis, labels_from_data = TRUE, pattern = NULL) {
    ## giftis is list of gifti images of length 2 (L hemi, R hemi)
    
    if (length(giftis) != 2) stop("giftis must have length 2")
    if (!(class(giftis[[1]]) == class(giftis[[2]]) && class(giftis[[1]]) == "gifti")) stop("elements of giftis must be of class gifti")
    if (isFALSE(all.equal(giftis[[1]]$data_info$Dim0, giftis[[2]]$data_info$Dim0))) stop("giftis must have same dimensions")
    if (!is.character(fold)) stop("fold must be class character")
    if (!is.logical(labels_from_data)) stop("labels_from_data must be logical")
    if (!(is.character(pattern) || is.null(pattern))) stop("pattern must be character or NULL")

    d <- lapply(giftis, function(x) abind::abind(x$data))  ## extract
    d <- abind::abind(d, along = 1)  ## L then R concat
    
    if (labels_from_data) {
        labs <- vapply(giftis[[1]]$data_meta, function(y) y[1, "vals"], character(1))
        hemis_are_mismatched <- !identical(labs, vapply(giftis[[2]]$data_meta, function(y) y[1, "vals"], character(1)))
        if (hemis_are_mismatched) stop("Regressor labels are mismatched across hemispheres.")
        dimnames(d) <- list(vertex = NULL, trial = labs)
        if (!is.null(pattern)) d <- d[, grep(pattern, labs)]
    }

    parcellated <- setNames(vector("list", length(rois)), names(rois))
    for (roi_i in seq_along(rois)) {
        which_parcels <- which(atlas$key$parcel %in% rois[[roi_i]])
        is_roi <- atlas$data %in% which_parcels
        d_roi <- d[is_roi, ]
        parcellated[[roi_i]] <- d_roi
    }
    
    ## get "bad list" of vertices
    good_vertices <- lapply(parcellated, function(x) apply(x, 1, var) > 0)    

    out <- list(list(data = parcellated, good_vertices = good_vertices, labels = labs))  ## nested
    names(out) <- fold

    ## merge lists element-wise into nested list
    #out <- c(Map(function(x, y) list(data = x, good_vertices = y), parcellated, good_list), labels = list(labs))

    class(out) <- c("parcellated_image", "list")

    out

}




## read atlas

## see here for more: https://github.com/mcfreund/psychomet/tree/master/in

read_atlas <- function(roiset = "schaefer", dir_atlas = "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/") {

    if (!roiset %in% c("schaefer", "mmp")) stop("roiset not configured")
    
    atlas <- list()

    if (roiset == "schaefer") {

        atlas$data <-
            c(
            gifti::read_gifti(file.path(dir_atlas, "Schaefer2018_400Parcels_7Networks_order_10K_L.label.gii"))$data[[1]],
            gifti::read_gifti(file.path(dir_atlas, "Schaefer2018_400Parcels_7Networks_order_10K_R.label.gii"))$data[[1]] + 200
            )
        atlas$key <- data.table::fread(here::here("in", "atlas-key_schaefer400-07.csv"))

    } else if (roiset == "mmp") {

    }

    atlas
    
}
