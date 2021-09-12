```{r}
library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel)
library(foreach)
library(mvnfast)
library(abind)
library(mikeutils)
library(lme4)
library(data.table)
library(ggplot2)
library(viridis)
library(colorspace)
library(grid)
library(gridExtra)
library(cowplot)

theme_set(theme_bw(base_size = 14))
theme_update(
  strip.background = element_rect(fill = "transparent", color = "transparent"),
)

vec <- function(x, diag = FALSE, ...) x[lower.tri(x, diag = FALSE), ...]


## build design matrices ----

labels <- rbind(
  expand.grid(
    target = c("black", "green", "pink", "yellow"),
    distractor = toupper(c("black", "green", "pink", "yellow")),
    stringsAsFactors = FALSE
  )
)
labels$congruency <- ifelse(labels$target != tolower(labels$distractor), "I", "C")
labels <- labels %>% arrange(target, distractor)


labels$count1 <- c(3, 0, 3, 3, 3, 6, 0, 0, 0, 3, 3, 3, 0, 3, 0, 6)
labels$count2 <- c(6, 3, 0, 0, 0, 3, 3, 3, 3, 0, 6, 0, 3, 0, 3, 3)

labels$run <- NA
labels$run[labels$count2 == 0] <- "run1"
labels$run[labels$count1 == 0] <- "run2"
labels$run[labels$count1 > 0 & labels$count2 > 0] <- "both"


trials <- 
  
  rbind(
    data.frame(
      target = rep(labels$target, labels$count1),
      distractor = rep(labels$distractor, labels$count1),
      run = 1
    ),
    
    data.frame(
      target = rep(labels$target, labels$count2),
      distractor = rep(labels$distractor, labels$count2),
      run = 2
      )
  )

trials$stimulus <- paste0(trials$target, "_", trials$distractor)
trials$congruency <- ifelse(trials$target != tolower(trials$distractor), "I", "C")
trials <- trials %>% group_by(run) %>% mutate(stimulus_trial = paste0(stimulus, "_", 1:n()))

# trials %>% filter(run == 1) %>% select(-stimulus) %>% table
# trials %>% filter(run == 1) %>% select(-stimulus, -distractor) %>% table
# trials %>% filter(run == 1)  %>% select(-stimulus, -target) %>% table
# 
# trials %>% 
#   filter(run == 1) %>%
#   filter(!(target == "pink" & distractor == "GREEN")) %>%
#   filter(!(target == "black" & distractor == "YELLOW")) %>%
#   select(-congruency, -stimulus) %>%
#   table
# 
# trials %>% 
#   filter(run == 2) %>%
#   filter(!(target == "green" & distractor == "PINK")) %>%
#   filter(!(target == "yellow" & distractor == "BLACK")) %>%
#   select(-congruency, -stimulus) %>%
#   table


## rsaX matrices (full matrix): across both runs (i.e., all stim):

X_s_all <- model.matrix(~ 0 + stimulus, trials)
colnames(X_s_all) <- gsub("stimulus", "", colnames(X_s_all))
rsX_all <- list(
  target = 
    tcrossprod(crossprod(X_s_all, model.matrix(~ 0 + target, trials)) > 0),
  distractor = 
    tcrossprod(crossprod(X_s_all, model.matrix(~ 0 + distractor, trials)) > 0),
  incongruency = tcrossprod(crossprod(X_s_all, cbind(0, model.matrix(~ 0 + congruency, trials)[, 2])) > 0),
  diag = tcrossprod(crossprod(X_s_all, X_s_all) > 0)
)




## X matrices: lists of length 2, one design matrix per run

X_s <- vector("list", 2)  ## stim
X_l <- vector("list", 2)  ## trial
X_t <- vector("list", 2)  ## target
X_d <- vector("list", 2)  ## distractor
X_i <- vector("list", 2)  ## incongruency

# rsx <- vector("list", 2)  ## rsa design matrix (off-diagonal)
rsX <- vector("list", 2)  ## rsa similarity matrices (nested list; full matrix)

for (run_i in 1:2) {
  # run_i = 1
  
  ## build design matrices:
  
  X_s[[run_i]] <- model.matrix(~ 0 + stimulus, trials %>% filter(run == run_i))
  X_l[[run_i]] <- model.matrix(~ 0 + stimulus_trial, trials %>% filter(run == run_i))
  X_t[[run_i]] <- model.matrix(~ 0 + target, trials %>% filter(run == run_i))
  X_d[[run_i]] <- model.matrix(~ 0 + distractor, trials %>% filter(run == run_i))
  X_i[[run_i]] <- model.matrix(~ 0 + congruency, trials %>% filter(run == run_i))
  X_i[[run_i]][, 1] <- 0  ## set congruent level to zero
  
  ## remove odd prefix from colname:
  
  colnames(X_s[[run_i]]) <- gsub("stimulus", "", colnames(X_s[[run_i]]))
  colnames(X_l[[run_i]]) <- gsub("stimulus_trial", "", colnames(X_l[[run_i]]))
  colnames(X_t[[run_i]]) <- gsub("target", "", colnames(X_t[[run_i]]))
  colnames(X_d[[run_i]]) <- gsub("distractor", "", colnames(X_d[[run_i]]))
  colnames(X_i[[run_i]]) <- gsub("congruency", "", colnames(X_i[[run_i]]))
  
  ## build RSA mats
  
  rsX[[run_i]] <- list(
    target = tcrossprod(crossprod(X_s[[1]], X_t[[1]])) > 0,
    distractor = tcrossprod(crossprod(X_s[[1]], X_d[[1]])) > 0,
    incongruency = tcrossprod(crossprod(X_s[[1]], X_i[[1]])) > 0
  )
  
}


rsx <- cbind(intercept = 1, 1 - do.call(cbind, lapply(rsX[[1]], vec)))





## plot ----



symmat4ggplot <- function(R, var.names = c("v1", "v2"), val.name = "value") {
  
  ## make factors for row and column labels
  dn <- dimnames(R)
  if (is.null(dn)) {
    dn <- setNames(list(paste0("cell_", 1:nrow(R)), paste0("cell_", 1:ncol(R))), var.names)
  } else {
    names(dn) <- var.names  
  }
  
  labels <- expand.grid(dn, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = TRUE)
  labels[[2]] <- factor(labels[[2]], levels = rev(levels(labels[[2]])))
  
  r <- c(R)
  
  cbind(labels, setNames(as.data.frame(r), val.name))
  
}



plotmat <- function(x, title, axis.text.size = 2, strip.text.size = 1.5, x_colors, y_colors, x_words, y_words) {
  
  # cols <- ifelse(as.character(labels$color) == "white", "grey50", as.character(labels$color))
  
  x %>%
    
    mutate(
      
      x = paste0(word1, "_", color1),
      y = paste0(word2, "_", color2),
      
      x = factor(x, levels = unique(x)),
      y = factor(y, levels = rev(unique(y)))
      
    ) %>%
    
    ggplot(aes(x, y, fill = value)) +
    geom_raster() +
    
    scale_fill_viridis_c(option = "magma") +
    scale_x_discrete(labels = x_words) +
    scale_y_discrete(labels = y_words) +
    
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_text(
        color = x_colors, face = "bold", size = rel(axis.text.size),
        angle = 90, hjust = 1, vjust = 0.5
      ),
      axis.text.y = element_text(color = y_colors, face = "bold", size = rel(axis.text.size)),
      axis.title = element_blank(),
      plot.title = element_text(size = rel(2)),
      panel.spacing = unit(0, "lines"), 
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = rel(strip.text.size))
    ) +
    
    labs(title = title)
  
  
}


labels$stimulus <- paste0(labels$target, "_", labels$distractor)

rsX_all$target_cross <- rsX_all$target[labels$stimulus[labels$count1 > 0], labels$stimulus[labels$count2 > 0]]
rsX_all$distractor_cross <- rsX_all$distractor[labels$stimulus[labels$count1 > 0], labels$stimulus[labels$count2 > 0]]
rsX_all$incongruency_cross <- rsX_all$incongruency[labels$stimulus[labels$count1 > 0], labels$stimulus[labels$count2 > 0]]

x_cross <- cbind(
  target = vec(rsX_all$target_cross),
  distractor = vec(rsX_all$distractor_cross),
  incongruency = vec(rsX_all$incongruency_cross)
)

x_full <- cbind(
  target = vec(rsX_all$target),
  distractor = vec(rsX_all$distractor),
  incongruency = vec(rsX_all$incongruency)
)


vif <- function(x) {
  vifs <- car::vif(lm(1:nrow(x) ~ ., as.data.frame(x)))
  data.frame(vifs, model = names(vifs))
}

vif(x_full)
vif(x_cross)

p_target <-
  rsX_all$target_cross %>%
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1")) %>%
  separate(v2, c("word2", "color2")) %>%
  plotmat(
    "target", 
    axis.text.size = 1, 
    x_colors = labels$target[labels$count1 > 0], 
    y_colors = rev(labels$target[labels$count2 > 0]),
    x_words = labels$distractor[labels$count1 > 0], 
    y_words = rev(labels$distractor[labels$count2 > 0])
  )

p_distractor <- 
  rsX_all$distractor_cross %>%
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1")) %>%
  separate(v2, c("word2", "color2")) %>%
  plotmat(
    "distractor", 
    axis.text.size = 1, 
    x_colors = labels$target[labels$count1 > 0], 
    y_colors = rev(labels$target[labels$count2 > 0]),
    x_words = labels$distractor[labels$count1 > 0], 
    y_words = rev(labels$distractor[labels$count2 > 0])
  )


p_incongruency <- 
  rsX_all$incongruency_cross %>%
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1")) %>%
  separate(v2, c("word2", "color2")) %>%
  plotmat(
    "incongruency", 
    axis.text.size = 1, 
    x_colors = labels$target[labels$count1 > 0], 
    y_colors = rev(labels$target[labels$count2 > 0]),
    x_words = labels$distractor[labels$count1 > 0], 
    y_words = rev(labels$distractor[labels$count2 > 0])
  )


p_pc50_cross <- arrangeGrob(
  p_target,
  p_distractor,
  p_incongruency,
  nrow = 1
  # layout_matrix = rbind(c(1:3, NA), 4:7)
)

ggsave(
  "C:/Users/mcf/Box/global/docs/papers/diss_proposal/figs/pc50_cross.pdf",
  p_pc50_cross,
  width = 4*3, height = 4, unit = "in", dev = "pdf"
)






## simulate ----

G <- 0*rsX_all$target + rsX_all$diag
# G[] <- 1
# diag(G) <- diag(G) + 1


get_coding_strengths <- function(
  covariance, 
  n_vertex = 100, 
  sd_B = 100, 
  sd_E = 1,
  .rsx = rsx_cv
) {
  
  n_trial1 <- nrow(X_s[[1]])
  n_trial2 <- nrow(X_s[[2]])
  stimuli1 <- colnames(X_s[[1]])
  stimuli2 <- colnames(X_s[[2]])
  n_stimuli <- nrow(G)
  
  ## generate data:
  
  B <- t(mvnfast::rmvn(n_vertex, mu = numeric(n_stimuli), sigma = G))
  E1 <- matrix(rnorm(n_trial1*n_vertex), nrow = n_trial1)
  E2 <- matrix(rnorm(n_trial2*n_vertex), nrow = n_trial2)
  
  Y1 <- X_s[[1]] %*% B[match(stimuli1, colnames(G)), ]*sd_B + E1*sd_E
  Y2 <- X_s[[2]] %*% B[match(stimuli2, colnames(G)), ]*sd_B + E2*sd_E
  
  ## estimate coefs:
  
  B1_hat <- coef(.lm.fit(X_s[[1]], Y1))
  B2_hat <- coef(.lm.fit(X_s[[2]], Y2))
  
  
  ## RSA:
  
  B_hat <- abind(t(B1_hat), t(B2_hat), along = 3)
  
  Y <- cvdist(B_hat, condition.names = paste0("tempname", 1:ncol(B_hat)))
  y <- vec(Y)
  
  image(cor(B_hat[, , 1], B_hat[, , 2]))
  
  b <- coef(.lm.fit(.rsx, y))
  names(b) <- colnames(.rsx)
  
  b
  
}

## how to build model RDM?
##  - need to check that rows/cols all line up
##  - issue with within-item (congr) vs btwn item (incon) comparison... model out?

rsx_t <- tcrossprod(
  crossprod(X_s[[1]], X_t[[1]]),
  crossprod(X_s[[2]], X_t[[2]])
) > 0
image(t(rsx_t))
rsx_d <- tcrossprod(
  crossprod(X_s[[1]], X_d[[1]]),
  crossprod(X_s[[2]], X_d[[2]])
) > 0

rsx_i <- tcrossprod(
  crossprod(X_s[[1]], X_i[[1]]),
  crossprod(X_s[[2]], X_i[[2]])
) > 0

rsx_stim <- tcrossprod(
  crossprod(X_s[[1]], X_i[[1]]),
  crossprod(X_s[[2]], X_i[[2]])
) > 0

X_s1_congr <- X_s[[1]]
X_s2_congr <- X_s[[2]]
image(X_s1_congr)
X_s1_congr[, !colnames(X_s[[1]]) %in% colnames(X_s[[2]])] <- 0
X_s2_congr[, !colnames(X_s[[2]]) %in% colnames(X_s[[1]])] <- 0

image(X_s2_congr)
rsx_stim <- tcrossprod(crossprod((1 - X_s1_congr), (1 - X_s2_congr)))
image(rsx_stim)



rsx_cv <- cbind(
  intercept = 1,
  target = 1 - vec(t(rsx_t)),
  distractor = 1 - vec(t(rsx_d)),
  incongruency = 1 - vec(t(rsx_i))
  )
# rsx_cv == rsx

image(t(rsX[[1]]$incongruency))
image(rsx_t)
image(rsx_d)
image(rsx_i)

```


