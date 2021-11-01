## system info

n_core <- parallel::detectCores()

## paths

dir_atlas <- "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/ATLASES/"


## image, design, analysis info

runs <- c("run1", "run2")
hemis <- c("L", "R")
sessions <- c("baseline", "proactive", "reactive")
sesss <- c("Bas", "Pro", "Rea")
waves <- c("wave1", "wave2")
wave_dir_image <- c(wave1 = "HCP_SUBJECTS_BACKUPS", wave2 = "DMCC_Phase3")
wave_dir_evts <- c(wave1 = "DMCC2", wave2 = "DMCC3")
colors_bias <- c("blue", "red", "purple", "white")
colors_pc50 <- c("black", "green", "pink", "yellow")
words_bias <- toupper(colors_bias)
words_pc50 <- toupper(colors_pc50)
stimuli_bias <- apply(expand.grid(colors_bias, words_bias), 1, paste0, collapse = "_")
stimuli_pc50 <- apply(expand.grid(colors_pc50, words_pc50), 1, paste0, collapse = "_")
trialtypes <- c(stimuli_bias, stimuli_pc50)


n_vertex <- 20484  ## surface hcp mesh
n_tr <- c(
  baseline  = 540,
  proactive = 540,
  reactive  = 590
)  ## number of tr per subj*run
n_trial <- c(
  baseline = 108,
  proactive = 108,
  reactive = 120
  )  ## number of trials (events) per subj*run
n_run <- 2
n_session <- 3
#n_trialtype <- length(trialtypes)  ## across both runs


## ROIs

core32 <- c(
  99, 127, 129, 130, 131, 132, 137, 140, 141, 142, 148, 163, 165, 182, 186, 300, 332, 333, 334, 335, 336, 337, 340, 345, 
  349, 350, 351, 352, 354, 361, 365, 387
)