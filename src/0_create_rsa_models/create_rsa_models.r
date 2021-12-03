library(here)
library(data.table)
library(magrittr)
source(here::here("src", "stroop-rsa-pc.R"))

## build label vectors

target <- gsub("[A-Z]", "", ttypes$all)
distractor <- gsub("[a-z]", "", ttypes$all)
incongruency <- as.character(target != tolower(distractor))

## build feature vectors (indicators)

ftrs_target <- indicator_matrix(target)
ftrs_distractor <- indicator_matrix(distractor)
ftrs_incongruency <- cbind(indicator_matrix(incongruency)[, "TRUE"])
ftrs_conjunction <- indicator_matrix(ttypes$all)
rownames(ftrs_target) <- ttypes$all
rownames(ftrs_distractor) <- ttypes$all
rownames(ftrs_incongruency) <- ttypes$all
rownames(ftrs_conjunction) <- ttypes$all

## make models for crcor: cross-run correlation matrices (with downsampling)

models_cor <- c("target", "distractor", "incongruency", "conjunction")
subsets <- c("all", "bias", "pc50")

l <- list(
    target = ftrs_target, 
    distractor = ftrs_distractor, 
    incongruency = ftrs_incongruency, 
    conjunction = ftrs_conjunction
    )

for (ses in sessions) {
    for (mod in models_cor) {
        for (sub in subsets) {

            conditions_run1 <- ttypes_by_run[[ses]]$run1
            conditions_run2 <- ttypes_by_run[[ses]]$run2

            if (sub == "bias") {
                conditions_run1 <- intersect(ttypes_by_run[[ses]]$run1, ttypes$bias)
                conditions_run2 <- intersect(ttypes_by_run[[ses]]$run2, ttypes$bias)
            } else if (sub == "pc50") {
                conditions_run1 <- intersect(ttypes_by_run[[ses]]$run1, ttypes$pc50)
                conditions_run2 <- intersect(ttypes_by_run[[ses]]$run2, ttypes$pc50)
            }

            x <- l[[mod]]
            X <- tcrossprod(x)[conditions_run1, conditions_run2]
            image(X)
            X <- as.data.table(X, keep.rownames = "run1")
            fwrite(X, here("out", "rsa_models", paste0("model_similarity_", mod, "_", ses, "_", sub, ".csv")))
            

        }
    }
}


## make models for cveuc: cross-validated euclidean (no downsampling)
models_dis <- c("target", "distractor", "incongruency")

for (ses in sessions) {
    for (mod in models_dis) {
        for (sub in subsets) {

            conditions_run1 <- ttypes_by_run[[ses]]$run1
            conditions_run2 <- ttypes_by_run[[ses]]$run2

            if (sub == "bias") {
                conditions_run1 <- intersect(conditions_run1, ttypes$bias)
                conditions_run2 <- intersect(conditions_run2, ttypes$bias)
            } else if (sub == "pc50") {
                conditions_run1 <- intersect(conditions_run1, ttypes$pc50)
                conditions_run2 <- intersect(conditions_run2, ttypes$pc50)
            }
            
            x <- l[[mod]]
            x1 <- x[conditions_run1, , drop = FALSE]
            x2 <- x[conditions_run2, , drop = FALSE]
            X <- cvdist(t(x1), t(x2))
            X <- X / abs(max(X)) ## scale so max val is 1
            image(X)
            dimnames(X) <- list(conditions_run1, conditions_run2)
            X <- as.data.table(X, keep.rownames = "run1")
            fwrite(X, here("out", "rsa_models", paste0("model_cvdistance_", mod, "_", ses, "_", sub, ".csv")))

        }
    }
}
