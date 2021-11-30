read_design <- function(subject, session, wave, glmname, run_i, signal_only) {

    dir_xmat <- here::here("out", "glms", subject, wave, "RESULTS/Stroop", paste0(session, "_", glmname))
    X <- mikeutils::read_xmat(file.path(dir_xmat, paste0("X.xmat_run", run_i ,".1D")))
    X <- X[rowSums(X) > .Machine$double.eps, ]  ## censor
    if (signal_only) {
        idx_baseline <- grep("Pol|movregs", colnames(X))
        X <- X[, -idx_baseline]  ## extract signal columns
    }
    X

}

estimate_distances <- function(B, measure, n_resamples, expected_min) {

    if (measure == "cveuc") {  ## cross-validated euclidean
        cvdist(
            average(B$run1, g = colnames(B$run1)),
            average(B$run2, g = colnames(B$run2)),
            scale = TRUE
        )
    } else if (measure == "crcor") {    ## cross-run correlation (with downsampling)
        crcor(B$run1, B$run2, n_resamples = n_resamples, expected_min = expected_min)
    }
}

regress_distances <- function(distances, design, measure, outcome_name = measure) {

    ## preproc:
    if (measure == "cveuc") {
        d <- vec(distances)
        Y <- d / sd(d)
    } else if (measure == "crcor") {
        d <- c(distances)
        Y <- atanh(d)  ## Fisher's z transform
    }

    Ahat <- coef(.lm.fit(x = design, y = Y))
    Ahat %>% tidy_model(terms = colnames(design), outcomes = outcome_name)

}


read_designs <- function(subjects, session, wave, glmname, signal_only) {
    X <- enlist(subjects)
    for (sub in subjects) {
        X[[sub]] <- enlist(runs)
        for (run_i in 1:2) {
            X[[sub]][[run_i]] <- read_design(sub, session, wave, glmname, run_i, signal_only = signal_only)
        }
    }    
    X
}


create_M <- function(ttypes) {

    p <- length(ttypes)
    m <- c("baseline", "diagonal", "distractor", "incongruency", "target")

    targets <- gsub("[A-Z]", "", ttypes)
    distractors <- gsub("[a-z]", "", ttypes)
    congruencies <- ifelse(targets == tolower(distractors), "congr", "incon")
    
    M_baseline <- matrix(1, nrow = p, ncol = p)
    M_diagonal <- diag(p)
    M_distractor <- tcrossprod(indicator_matrix(distractors))
    M_incongruency <- tcrossprod(indicator_matrix(congruencies)[, "incon"])
    M_target <- tcrossprod(indicator_matrix(targets))

    cbind(
        baseline = c(M_baseline), 
        diagonal = c(M_diagonal), 
        distractor = c(M_distractor), 
        incongruency = c(M_incongruency), 
        target = c(M_target)
        )
    
}



simulate_experiments <- function(
    X_list, Z_list, Q_sim, Q_dis,
    .glmname = glmname, .ses = ses, .wav = wav, .ttype_subset = ttype_subset, .n_experiments = n_experiments, 
    .mu = mu, .sigma = sigma, .v = v, .r = r, .n_resamples = n_resamples, 
    .expected_min = expected_min[paste0(.ses, "_", .ttype_subset)], est_crcor = TRUE, .n_cores = n_cores
    ) {

    stopifnot(identical(names(X_list), names(Z_list)))
    stopifnot(!is.null(names(X_list)))
    subjects <- names(X_list)
    
    ## trial labels

    tinfo <- read_trialinfo()[subj %in% subjects & wave %in% .wav & session %in% .ses]
    setkey(tinfo, subj, run)  ## for quick subsetting within loop below

    ## simulation loop:
    
    set.seed(0)
    B_array <- replicate(length(subjects), t(MASS::mvrnorm(n = .v, mu = .mu, Sigma = .sigma)))  ## create betas  
    
    cl <- makeCluster(.n_cores, type = "FORK")        
    registerDoParallel(cl)
    res <- 
        foreach(subj_i = seq_along(subjects), .final = function(x) setNames(x, subjects), .inorder = TRUE) %:% 
        foreach(experiment_i = seq_len(.n_experiments), .inorder = FALSE) %dopar% {

            set.seed((subj_i-1)*.n_experiments + experiment_i)  ## uniqe seed per subject and experiment

            sub <- subjects[subj_i]
            X <- X_list[[subj_i]]
            Z <- Z_list[[subj_i]]
            B <- B_array[, , subj_i]

            Bhat <- enlist(runs)
            for (run_i in 1:2) {

                ## generate timeseries Y
                
                t <- nrow(X[[run_i]])  ## n tr (after censoring)
                E <- matrix(rnorm(.v*t, sd = sqrt(r)), ncol = .v)
                Y <- X[[run_i]] %*% B + E

                ## estimate coefficients Bhat
                
                idx_signal <- -grep("Pol|movregs", colnames(Z[[run_i]]))
                Bhat_i <- t(coef(.lm.fit(Z[[run_i]], Y))[idx_signal, ])

                tlabels <- tinfo[.(sub, run_i), item]  ## get trialtypes (add as colnames)
                if (glmname == "condition_1rpm") tlabels <- c(ttypes$bias, ttypes$pc50)
                colnames(Bhat_i) <- tlabels

                ttypes_to_get <- intersect(ttypes_by_run[[.ses]][[run_i]], ttypes[[.ttype_subset]])  ## extract ttypes
                Bhat_i <- Bhat_i[, colnames(Bhat_i) %in% ttypes_to_get]

                Bhat[[run_i]] <- Bhat_i

            }
            
            ahat <-
                Bhat %>%
                estimate_distances("cveuc") %>%
                regress_distances(Q_dis, "cveuc") %>%
                melt(id.vars = "term")

            if (est_crcor) {
                ahat_sim <-
                    Bhat %>%
                    estimate_distances("crcor", n_resamples = n_resamples, expected_min = .expected_min) %>%
                    regress_distances(Q_sim, "crcor") %>%
                    melt(id.vars = "term")
                ahat <- rbind(melt(ahat_sim, id.vars = "term"), ahat)
            }

            ahat

    }
    stopImplicitCluster()

    res_unnested <- lapply(res, rbindlist, id = "experiment")
    
    rbindlist(res_unnested, id = "subj")

}

