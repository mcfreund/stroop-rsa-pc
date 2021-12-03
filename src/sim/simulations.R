library(colorout)
library(here)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel)
library(foreach)

source(here("src", "stroop-rsa-pc.R"))
source(here("src", "sim", "functions.R"))

theme_set(theme_bw(base_size = 10))



## constants ----

subjlist <- "all_retest"
subjects <- fread(here("out", paste0("subjlist_", subjlist, ".txt")))[[1]]
glmname <- "lsall_1rpm"
wav <- "wave1"
n_experiments <- 100  ## number simulations per subject
n_cores <- 16
p <- length(ttypes$all)
Mvec <- create_M(ttypes$all)  ## components of true similarity structure
mu <- rep(0, p)  ## true mean activations
v <- 100  ## number vertices
n_resamples <- 1E2  ## for crcor
r <- 10^-c(4:-1)  ## SNRs


## read xmats ----

## read trial orders (for lsall):

tinfo <- read_trialinfo()[subj %in% subjects & wave %in% wav, c("subj", "session", "run", "trial_num", "item")]
setkey(tinfo, subj, session, run, trial_num)  ## for quick subsetting


xmat_condition_master <- enlist(sessions)
xmat_lsall_master <- enlist(sessions)
xmat_simil_master <- enlist(sessions)
xmat_dist_master <- enlist(sessions)

for (ses in sessions) {

    ## load timeseries xmats
    xmat_condition_master[[ses]] <- read_designs(subjects, ses, wav, "condition_1rpm")
    xmat_lsall_master[[ses]] <- read_designs(subjects, ses, wav, "lsall_1rpm")
    
    ## rename "signal" (stimulus) cols:
    for (subject in subjects) {
        for (run_i in 1:2) {
            
            nms_lsall <- colnames(xmat_lsall_master[[ses]][[subject]][[run_i]])
            nms_lsall[grep("alltrials", nms_lsall)] <- tinfo[.(subject, ses, run_i), item]
            colnames(xmat_lsall_master[[ses]][[subject]][[run_i]]) <- nms_lsall

            nms_condition <- colnames(xmat_condition_master[[ses]][[subject]][[run_i]])
            nms_condition <- gsub("#0", "", nms_condition)
            colnames(xmat_condition_master[[ses]][[subject]][[run_i]]) <- nms_condition

        }
    }
    
    ## load RSA xmats
    xmat_simil_master[[ses]] <- enlist(names(ttypes))
    xmat_dist_master[[ses]] <- enlist(names(ttypes))

    for (ttype_subset in names(ttypes)) {
        xmat_simil_master[[ses]][[ttype_subset]] <- read_model_xmat("crcor", ses, ttype_subset)
        xmat_dist_master[[ses]][[ttype_subset]] <- read_model_xmat("cveuc", ses, ttype_subset)
    }

}




## false positive simulations ----

a <- c(1, 0, 0, 0, 0)  ## true model weights
sigma <- matrix(Mvec %*% a, ncol = p, nrow = p, dimnames = list(ttypes$all, ttypes$all))  ## true geometry
params <- expand.grid(
    ttype_subset = c("pc50", "bias"), #names(ttypes),
    ses = sessions,
    r = r,
    stringsAsFactors = FALSE
)

beg <- Sys.time()
dat <- vector("list", nrow(params))
for (param_i in seq_len(nrow(params))) {
    # param_i = 10
    
    ttype_subset <- params$ttype_subset[param_i]
    ses <- params$ses[param_i]

    result <- simulate_experiments(
        .X_list = xmat_condition_master[[ses]], 
        .Z_list = xmat_lsall_master[[ses]], 
        .Q_dis  = xmat_dist_master[[ses]][[ttype_subset]],
        .Q_sim  = xmat_simil_master[[ses]][[ttype_subset]], 
        .ttype_subset = ttype_subset,
        .r = params$r[param_i]
    )

    ## save
    
    result$session <- ses
    result$ttype_subset <- ttype_subset
    result$r <- params$r[param_i]

    dat[[param_i]] <- result
    
    print(param_i)

}
false_positive <- rbindlist(dat)
end <- Sys.time()
(end - beg)

saveRDS(false_positive, here("out", "sim", paste0("false_positive_", wav, ".RDS")))



## true positive simulations ----

a_list <- list(
    distractor_target                          = c(0.5, 0, 0.5/4, 0, 0.5/4),
    incongruency                               = c(0.5, 0, 0, 0.5/4, 0),
    distractor_target_incongruency             = c(0.5, 0, 0.5/4, 0.5/4, 0.5/4),
    distractor_target_incongruency_conjunction = c(0.5, 0.5/4, 0.5/4, 0.5/4, 0.5/4)
)
for (effect_i in seq_along(a_list)) {

    a <- a_list[[effect_i]]
    sigma <- matrix(Mvec %*% a, ncol = p, nrow = p, dimnames = list(ttypes$all, ttypes$all))  ## true geometry
    params <- expand.grid(
        ttype_subset = c("pc50", "bias"),
        ses = sessions,
        r = r,
        stringsAsFactors = FALSE
    )

    beg <- Sys.time()
    dat <- vector("list", nrow(params))
    for (param_i in seq_len(nrow(params))) {
        # param_i = 10
        
        ttype_subset <- params$ttype_subset[param_i]
        ses <- params$ses[param_i]

        result <- simulate_experiments(
            .X_list = xmat_condition_master[[ses]], 
            .Z_list = xmat_lsall_master[[ses]], 
            .Q_dis  = xmat_dist_master[[ses]][[ttype_subset]],
            .Q_sim  = xmat_simil_master[[ses]][[ttype_subset]], 
            .ttype_subset = ttype_subset,
            .r = params$r[param_i]
        )

        ## save
        
        result$session <- ses
        result$ttype_subset <- ttype_subset
        result$r <- params$r[param_i]

        dat[[param_i]] <- result
        rm(result)
        gc()
        
        print(paste0(param_i, "/", nrow(params)))

    }
    end <- Sys.time()
    (end - beg)

    saveRDS(rbindlist(dat), here("out", "sim", paste0("true_positive_", names(a_list[effect_i]), "_", wav, ".RDS")))
    rm(dat)
    gc()

}
