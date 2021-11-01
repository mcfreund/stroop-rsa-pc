## 
## defines functions for use in this project

## class def

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
