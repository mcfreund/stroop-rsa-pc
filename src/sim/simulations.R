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
n_experiments <- 1000  ## number simulations per subject
n_cores <- 26
p <- length(ttypes$all)
Mvec <- create_M(ttypes$all)  ## components of true similarity structure
a <- setNames(rep(0, ncol(Mvec)), colnames(Mvec))  ## true model weights
mu <- rep(0, p)  ## true mean activations
v <- 100  ## number vertices
n_resamples <- 1E3  ## for crcor
r <- seq(1E-6, 1, 1/4)  ## SNRs
 
## read all design matrices into nested list:

X_master <- enlist(sessions)
Z_master <- enlist(sessions)
Q_sim_master <- enlist(sessions)
Q_dis_master <- enlist(sessions)
for (session in sessions) {
    X_master[[session]] <- read_designs(subjects, session, wav, "condition_1rpm", signal_only = TRUE)
    Z_master[[session]] <- read_designs(subjects, session, wav, glmname, signal_only = FALSE)
    Q_sim_master[[session]] <- enlist(names(ttypes))
    Q_dis_master[[session]] <- enlist(names(ttypes))
    for (ttype_subset in names(ttypes)) {
        Q_sim_master[[session]][[ttype_subset]] <- read_model_xmat("crcor", session, ttype_subset)
        Q_dis_master[[session]][[ttype_subset]] <- read_model_xmat("cveuc", session, ttype_subset)
    }
}



## false positive simulations ----

a_baseline <- seq(0, 0.99, 0.9/4)
params <- expand.grid(
    ttype_subset = "bias", #c("bias", "pc50", "all"),
    ses = "reactive",
    r = r,
    a_baseline = a_baseline,
    stringsAsFactors = FALSE
)
params$a_diagonal <- 1 - a_baseline




beg <- Sys.time()

dat <- vector("list", nrow(params))
for (param_i in seq_len(nrow(params))) {
    # param_i = 1
    
    ttype_subset <- params$ttype_subset[param_i]
    ses <- params$ses[param_i]

    ##  population parameters

    .r <- params$r[param_i]
    .a <- c(params$a_baseline[param_i], params$a_diagonal[param_i], 0, 0, 0)    
    .sigma <- matrix(Mvec %*% .a, ncol = p, nrow = p, dimnames = list(ttypes$all, ttypes$all))  ## true geometry

    ## simulate
    
    result <- simulate_experiments(
        X_list = X_master[[ses]], 
        Z_list = Z_master[[ses]], 
        Q_sim = Q_sim_master[[ses]][[ttype_subset]], 
        Q_dis = Q_dis_master[[ses]][[ttype_subset]],
        .ses = ses,
        .ttype_subset = ttype_subset,
        .sigma = .sigma,
        .r = .r
    )
    
    ## save
    
    result$session <- ses
    result$ttype_subset <- ttype_subset
    result$a_baseline <- params$a_baseline[param_i]
    result$r <- .r

    dat[[param_i]] <- result

}
false_positive <- rbindlist(dat)

end <- Sys.time()
