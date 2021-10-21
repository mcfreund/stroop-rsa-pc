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

theme_set(theme_bw(base_size = 14))
theme_update(
  strip.background = element_rect(fill = "transparent", color = "transparent"),
)


## build 'label' dataframes ----
## that hold experimental condition information

labels <- rbind(
  # expand.grid(
  #   target = c("blue", "purple", "red", "white"),
  #   distractor = toupper(c("blue", "purple", "red", "white"))
  # ),
  expand.grid(
    target = c("black", "green", "pink", "yellow"),
    distractor = toupper(c("black", "green", "pink", "yellow"))
  )
)
labels$congruency <- ifelse(labels$target != tolower(labels$distractor), "I", "C")
labels <- lapply(labels, as.factor)
labels_a <- data.frame(
  target = paste0(labels$target, "_", labels$congruency),
  distractor = paste0(labels$distractor, "_", labels$congruency),
  stringsAsFactors = TRUE
  )

## transform labels to dummy codes:

is_level <- model.matrix(~ . + 0, labels, contrasts.arg = lapply(labels, contrasts, contrasts = FALSE)) == 1
is_level_a <- model.matrix(~ . + 0, labels_a, contrasts.arg = lapply(labels_a, contrasts, contrasts = FALSE)) == 1



## functions ----
## NB: no downsampling yet, just aggregation

vec <- function(x, diag = FALSE, ...) x[lower.tri(x, diag = FALSE), ...]
offdiag <- function(x, ...) x[col(x) != row(x), ...]
lapply_names <- function(nms, fun, ...) setNames(lapply(nms, fun, ...), nms)



get_X <- function(n_reps_per_cond, conds) {
    ## generate timeseries design matrices X
    
    # n_reps_per_cond <- 10
    # conds <- paste0(labels$target, "_", labels$distractor)
    
    ## get X
    
    cond_order <- rep(conds, n_reps_per_cond)
    X <- model.matrix(~ 0 + cond_order)
    colnames(X) <- gsub("cond_order", "", colnames(X))
    X <- X[, conds]  ## rearrange
    
    ## get X_trial
    
    cond_order_trial <- paste0(cond_order, "_", 1:length(cond_order))
    X_trial <- model.matrix(~ 0 + cond_order_trial)
    colnames(X_trial) <- gsub("^.*cond_order_trial", "", colnames(X_trial))
    X_trial <- X_trial[, cond_order_trial]  ## rearrange!
    
    ## get X_models, bases for models
    
    X_distractor <- X %*% is_level[, grep("distractor", colnames(is_level))]
    # X_incongruency <- X %*% is_level[, grep("congruency", colnames(is_level))]
    X_incongruency <- X %*% cbind(FALSE, is_level[, grep("congruencyI", colnames(is_level))])
    X_target <- X %*% is_level[, grep("target", colnames(is_level))]
    
    
    list(
      stimuli = X, 
      trial = X_trial, 
      distractor = X_distractor,
      target = X_target,
      incongruency = X_incongruency
    )
    
}



get_A <- function(X_stimuli) {

  A_stimuli <- sweep(X_stimuli, 2, colSums(X_stimuli), "/")
  
  A_target <- is_level_a[, grep("target", colnames(is_level_a))]
  colnames(A_target) <- gsub("target", "", colnames(A_target))
  A_target <- sweep(A_target, 2, colSums(A_target), "/")
  
  A_distractor <- is_level_a[, grep("distractor", colnames(is_level_a))]
  colnames(A_distractor) <- gsub("distractor", "", colnames(A_distractor))
  A_distractor <- sweep(A_distractor, 2, colSums(A_distractor), "/")
  
  A_target <- A_stimuli %*% A_target
  A_distractor <- A_stimuli %*% A_distractor
  
  list(
    stimuli = A_stimuli, 
    target = A_target, 
    distractor = A_distractor
  )
  
}




get_rsa_x <- function(name_A, as_xmat = TRUE, all_offdiag = FALSE) {
  ## get RSA regression design matrices
  
  rsa_x_bar_t <- tcrossprod(crossprod(A[[name_A]], X$target))
  rsa_x_bar_d <- tcrossprod(crossprod(A[[name_A]], X$distractor))
  rsa_x_bar_i <- tcrossprod(crossprod(A[[name_A]], X$incongruency))
  
  
  if (as_xmat) {
    
    if (all_offdiag == FALSE) {
      
      cbind(
        intercept = 1,
        target = vec(rsa_x_bar_t),
        distractor = vec(rsa_x_bar_d),
        incongruency = vec(rsa_x_bar_i)
      )
      
    } else {
      
      cbind(
        intercept = 1,
        target = offdiag(rsa_x_bar_t),
        distractor = offdiag(rsa_x_bar_d),
        incongruency = offdiag(rsa_x_bar_i)
      )
      
    }
    
  } else {
    
    list(
      target = rsa_x_bar_t,
      distractor = rsa_x_bar_d,
      incongruency = rsa_x_bar_i
    )
    
  }
  
}


get_coding_strengths <- function(
  name, S, 
  .n_vertex =  n_vertex, 
  models = c("intercept", "target", "distractor", "incongruency"), 
  rsatype = "cr",
  sd_B = 1, 
  sd_E = 1
  ) {
  
  n_trial <- nrow(X$stimuli)
  n_stimuli <- ncol(X$stimuli)
  
  ## generate data:
  
  # B <- t(mvnfast::rmvn(.n_vertex, mu = numeric(n_stimuli), sigma = S))
  B <- t(matrix(rnorm(.n_vertex*n_stimuli), ncol = n_stimuli) %*% S)
  E1 <- matrix(rnorm(n_trial*.n_vertex), nrow = n_trial)
  E2 <- matrix(rnorm(n_trial*.n_vertex), nrow = n_trial)
  
  Y1 <- X$stimuli %*% B*sd_B + E1*sd_E
  Y2 <- X$stimuli %*% B*sd_B + E2*sd_E
  
  ## estimate coefs:
  
  B1_hat <- coef(.lm.fit(X$trial, Y1))
  B2_hat <- coef(.lm.fit(X$trial, Y2))
  
  
  B1_hat_bar <- crossprod(A[[name]], B1_hat)  ## aggregate over trials (and possibly levels)
  B2_hat_bar <- crossprod(A[[name]], B2_hat)  ## aggregate over trials (and possibly levels)
  
  
  ## RSA:
  
  .rsa_x <- rsa_x[[name]][, models]
  
  if (rsatype == "cv") {
    
    B_hat <- abind(t(B1_hat_bar), t(B2_hat_bar), along = 3)
    B_hat <- -B_hat
    y <- vec(cvdist(B_hat, condition.names = colnames(B_hat)))

  } else if (rsatype == "cr") {
    
    y <- offdiag(atanh(cor(t(B1_hat_bar), t(B2_hat_bar))))
    
  }
  
  
  b <- coef(.lm.fit(.rsa_x, y))
  names(b) <- colnames(.rsa_x)
  
  b
  
}


## use ----

## get design and sim info:

n_core <- parallel::detectCores()
n_sim <- 1000
n_vertex <- 100
n_reps_per_cond <- 10
conds <- paste0(labels$target, "_", labels$distractor)
X <- get_X(n_reps_per_cond, conds)
A <- get_A(X$stimuli)
rsa_x <- lapply_names(c("stimuli", "target", "distractor"), get_rsa_x, all_offdiag = TRUE)
rsa_X_stimuli <- get_rsa_x("stimuli", as_xmat = FALSE)

snr <- c(0, 10^(-2:0))
# b <- list(
#   t = c(t = 1, d = 0, i = 0, diag = 1),
#   d = c(t = 0, d = 1, i = 0, diag = 1),
#   i = c(t = 0, d = 0, i = 1, diag = 1),
#   td = c(t = 1, d = 1, i = 0, diag = 1),
#   ti = c(t = 1, d = 0, i = 1, diag = 1),
#   di = c(t = 0, d = 1, i = 1, diag = 1),
#   tdi = c(t = 1, d = 1, i = 1, diag = 1)
# )


b <- list(
  c(t = 0.2, d = 0, i = 0, diag = 1),
  c(t = 0.2, d = 0.2, i = 0, diag = 1),
  c(t = 0.2, d = 0.4, i = 0, diag = 1),
  c(t = 0.2, d = 0.8, i = 0, diag = 1)
)

for (ii in seq_along(rsa_X_stimuli)) diag(rsa_X_stimuli[[ii]]) <- 0  ## remove diagonal

## check positive definitness:
for (geom_i in seq_along(b)) {
  b_i <- b[[geom_i]]
  S <-
    b_i["t"] * rsa_X_stimuli$target +
    b_i["d"] * rsa_X_stimuli$distractor +
    b_i["i"] * rsa_X_stimuli$incongruency +
    b_i["diag"] * diag(nrow(rsa_X_stimuli$target))
  chol(S)
}

## simulate:

l <- setNames(vector("list", length(names(b))), names(b))
for (geom_i in seq_along(b)) {
  
  
  b_i <- b[[geom_i]]
  S <- 
    b_i["t"] * rsa_X_stimuli$target +
    b_i["d"] * rsa_X_stimuli$distractor +   
    b_i["i"] * rsa_X_stimuli$incongruency + 
    b_i["diag"] * diag(nrow(rsa_X_stimuli$target))
  

  set.seed(0)
  
  cl <- makeCluster(n_core - 1)
  registerDoParallel(cl)
  
  res <-
    foreach(ii = snr, .inorder = FALSE, .combine = "rbind", .packages = c("data.table", "abind", "mikeutils")) %:%
    foreach(jj = seq_len(n_sim), .inorder = FALSE, .combine = "rbind") %dopar% {
      
      b_stimuli <- get_coding_strengths("stimuli", S, sd_B = ii)
      b_agg_t <- get_coding_strengths("target", S, models = c("intercept", "target", "incongruency"), sd_B = ii)
      b_agg_d <- get_coding_strengths("distractor", S, models = c("intercept", "distractor", "incongruency"), sd_B = ii)
      b_agg_t_cov <- get_coding_strengths("target", S, sd_B = ii)
      b_agg_d_cov <- get_coding_strengths("distractor", S, sd_B = ii)
      
      v <- 
        unlist(
          list(
            stimuli = b_stimuli, 
            agg_t_cov = b_agg_t_cov, 
            agg_d_cov = b_agg_d_cov, 
            agg_t = b_agg_t, 
            agg_d = b_agg_d
            )
          )
      
      data.table(b = v, variable = names(v), snr = ii)
      
      
    }
    
  stopCluster(cl)
  
  l[[geom_i]] <- res %>% separate(variable, into = c("modeltype", "term"), "\\.")
  
  
}



d <- rbindlist(l, idcol = "geometry")
# d$geometry <- factor(d$geometry, levels = names(b))
d$geometry <- factor(d$geometry)

# saveRDS(d, "C:/Users/mcf/Box/global/docs/papers/diss_proposal/misc/sim_results.RDS")




## plot:

mean_plusminus_sd <- function(x) {
  mu <- mean(x)
  sigma <- sd(x)
  data.frame(y = mu, ymax = mu + sigma, ymin = mu - sigma)
}


is_odd_comb <- 
  (d$term == "distractor" & d$modeltype == "agg_t_cov") |
  (d$term == "target" & d$modeltype == "agg_d_cov")


d_subset <- d[!is_odd_comb] 


d_subset %>%
  
  ggplot(aes(as.factor(snr), b, color = modeltype)) +
  
  stat_summary(
    fun.data = mean_plusminus_sd, geom = "pointrange", 
    position = position_dodge(width = 0.5),
    size = 1
    ) +

  facet_grid(rows = vars(geometry), cols = vars(term)) +
  scale_color_brewer(type = "qual", palette = 2) +
  theme(legend.position = c(0.6, 0.8)) +
  
  labs(x = "SNR", y = "mean b +/- SD")

d_subset %>%
  
  group_by(modeltype, term, geometry, snr) %>%
  mutate(
    sim = 1:n(),
    experiment = ntile(sim, 25)
    ) %>%
  group_by(term, modeltype, snr, geometry, experiment) %>%
  summarize(stat = t.test(b)$statistic) %>%
  summarize(stat = mean(stat)) %>%
  
  ggplot(aes(as.factor(snr), stat, color = modeltype)) +
  geom_hline(yintercept = 0) + 
  geom_line(aes(group = modeltype), size = 1.5) +
  
  facet_grid(rows = vars(geometry), cols = vars(term)) +
  scale_color_brewer(type = "qual", palette = 2) +
  theme(legend.position = c(0.6, 0.8)) +
  
  labs(x = "SNR", y = "mean t (N = 25)")


d_subset %>%
  
  group_by(modeltype, term, geometry, snr) %>%
  mutate(
    sim = 1:n(),
    experiment = ntile(sim, 25)
  ) %>% 
  
  mutate(modeltype = gsub("_d|_t", "", modeltype))

  group_by(term, snr, geometry, modeltype, experiment) %>%
  summarize(stat = t.test(b)$statistic) %>%
  summarize(stat = mean(stat)) %>%
  
  ggplot(aes(as.factor(snr), stat, color = modeltype)) +
  geom_hline(yintercept = 0) + 
  geom_line(aes(group = modeltype), size = 1.5) +
  
  facet_grid(rows = vars(geometry), cols = vars(term)) +
  scale_color_brewer(type = "qual", palette = 2) +
  theme(legend.position = c(0.6, 0.8)) +
  
  labs(x = "SNR", y = "mean t (N = 25)")

  
## for examining influence of distractor on target:  

d_subset %>%
  
  group_by(modeltype, term, geometry, snr) %>%
  mutate(
    sim = 1:n(),
    experiment = ntile(sim, 25)
  ) %>%
  group_by(term, modeltype, snr, geometry, experiment) %>%
  summarize(stat = t.test(b)$statistic) %>%
  summarize(stat = mean(stat)) %>%
  
  ggplot(aes(as.factor(snr), stat, color = geometry)) +
  geom_hline(yintercept = 0) + 
  geom_line(aes(group = geometry), size = 1.5) +
  
  facet_grid(rows = vars(modeltype), cols = vars(term)) +
  scale_color_viridis_d() +
  theme(legend.position = c(0.6, 0.8)) +
  
  labs(x = "SNR", y = "mean t (N = 25)")

d_subset %>%
  
  group_by(term, modeltype, snr, geometry) %>%
  summarize(stat = mean(b)) %>%
  
  ggplot(aes(as.factor(snr), stat, color = geometry)) +
  geom_hline(yintercept = 0) + 
  geom_line(aes(group = geometry), size = 1.5) +
  
  facet_grid(rows = vars(modeltype), cols = vars(term)) +
  scale_color_viridis_d() +
  theme(legend.position = c(0.6, 0.8)) +
  
  labs(x = "SNR", y = "mean b (N = 25)")





## are aggregation coefficients biased?
## does adding a covariate help?
## should models be fit at stimulus level anyway?

## for true incongruent coding
##    - aggregation is more powerful/(biased positive?) than stimulus-level, other ests not impacted
##    - no impact of covariate on aggregation
## for true target coding
##    - largest target t-stat in stimulus-level (esp good SNR), agg-t and agg-t-cov next, then agg-d
##      - perhaps a positive bias in mean beta in agg-t/agg-t-cov, but larger variance
##    - in incongruent estimates, positive bias in agg-d and agg-d-cov (no impact of cov) 
##    - in distractor estimates, negative bias in agg-d, but not agg-d-cov (or any others)
##    - in intercept, positive bias in agg-d, slight negative bias in agg-d-cov
## for true distractor coding
##    - symmetric with true target coding
## general notes:
##    - aggregation method is assisted (made less biased) by including covariate, as negaitve effect of target 
##      (or distractor) on other is reduced.
##    - maybe averaging incongruency across models will reduce bias... but still positive, since both target and
##      distractor INCREASE incongruency coding.
##    - with congruency model (not incongruency), different behavior between agg and agg cov incongruency estimate when 
##      simultaneous target and distractor coding? positive bias in agg cov, negative bias in agg?


## 


# res_d <- res
# res_d <- res_d %>%
#   
#   group_by(modeltype, term, e_var) %>%
#   mutate(
#     sim = 1:n(), 
#     experiment = ntile(sim, 25)
#     )
# 
# res_d %>%
# 
#   group_by(term, modeltype, e_var, experiment) %>%
# 
#   summarize(stat = t.test(b)$statistic, p = t.test(b)$p.value) %>%
#   summarize(hitrate = sum(p < 0.05)/n()) %>%
# 
#   ggplot(aes(as.factor(e_var), hitrate, color = modeltype)) +
#   geom_line(aes(group = modeltype), size = 1) +
#   facet_grid(cols = vars(term)) +
#   scale_color_brewer(type = "qual", palette = 2) +
#   scale_x_discrete(labels = round(unique(res$e_var), 0)) +
#   theme(legend.position = c(0.9, 0.8), axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
#   labs(x = "error variance", y = "hit rate")
# 




## assessing permutation test ----



n_perm <- 100
n_sim <- 100

# DEFINE PERM MATRIX WITH ROWS THAT GIVE SHUFFLES
# perms
set.seed(0)

perms <- data.frame(orig = colnames(X$stimuli))
setDT(perms)
perms <- perms %>% separate(orig, c("target", "distractor"), "_", remove = FALSE)
perms$congruency <- ifelse(perms$target == tolower(perms$distractor), "C", "I")
l <- vector("list", n_perm)
for (ii in seq_len(n_perm)) {
  
  l[[ii]] <- perms[, perm := sample(orig), distractor][, c("orig", "perm")]
  
} 



cl <- makeCluster(n_core - 1)
registerDoParallel(cl)

res <-
  foreach(ii = seq_len(n_perm), .inorder = FALSE, .combine = "rbind", .packages = c("data.table", "abind", "mikeutils")) %:%
  foreach(jj = seq_len(n_sim), .inorder = FALSE, .combine = "rbind") %dopar% {
    
    # perm_i <- l[[ii]]$perm
    # S_star <- S[perm_i, perm_i]
    
    S_star <- S
    
    b_agg_t <- get_coding_strengths(
      e_var = 50, "target", S = S_star, n_vertex = n_vertex, models = c("intercept", "target", "incongruency"),
      all_offdiag = FALSE
    )
    b_agg_t_cov <- get_coding_strengths(e_var = 50, "target", S = S_star, n_vertex = n_vertex)
    
    
    v <- unlist(list(agg_t_cov = b_agg_t_cov, agg_t = b_agg_t))
    
    data.table(b = v, variable = names(v))
    
    
  }

stopCluster(cl)



res <- res %>% separate(variable, into = c("modeltype", "term"), "\\.")

res %>%
  
  group_by(term, modeltype) %>%
  
  # summarize(stat = t.test(b)$statistic) %>%
  summarize(stat = mean(b)) %>%
  
  ggplot(aes(modeltype, stat, fill = modeltype)) +
  geom_col(position = position_dodge(width = 1)) +
  facet_grid(cols = vars(term)) +
  scale_fill_brewer(type = "qual", palette = 2) +
  theme(legend.position = c(0.9, 0.8), axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))


## getting unbiased stimulus set for i-c contrast ----

labels <- rbind(
  expand.grid(
    target = c("black", "green", "pink", "yellow"),
    distractor = toupper(c("black", "green", "pink", "yellow")),
    stringsAsFactors = FALSE
  )
)
labels$congruency <- ifelse(labels$target != tolower(labels$distractor), "I", "C")
labels <- labels %>% arrange(target, distractor)

# labels$run <- c(
#   "both", "run1", "run2", "run2",  ## BLACK word: alphabetical by color
#   "run2", "both", "run1", "run1",  ## GREEN
#   "run1", "run2", "both", "run2",  ## PINK
#   "run1", "run2", "run1", "both"   ## YELLOW
#   )

labels$count1 <- c(3, 0, 3, 3, 3, 6, 0, 0, 0, 3, 3, 3, 0, 3, 0, 6)
labels$count2 <- c(6, 3, 0, 0, 0, 3, 3, 3, 3, 0, 6, 0, 3, 0, 3, 3)

labels$run <- NA
labels$run[labels$count2 == 0] <- "run1"
labels$run[labels$count1 == 0] <- "run2"
labels$run[labels$count1 > 0 & labels$count2 > 0] <- "both"

trials1 <- data.frame(
  target = rep(labels$target, labels$count1),
  distractor = rep(labels$distractor, labels$count1)
)
trials1$congruency <- ifelse(trials1$target != tolower(trials1$distractor), "I", "C")


trials2 <- data.frame(
  target = rep(labels$target, labels$count2),
  distractor = rep(labels$distractor, labels$count2)
)
trials2$congruency <- ifelse(trials2$target != tolower(trials2$distractor), "I", "C")

table(trials1)
table(trials1 %>% select(-distractor))
table(trials1 %>% select(-target))


trials1 %>% 
  filter(!(target == "pink" & distractor == "GREEN")) %>%
  filter(!(target == "black" & distractor == "YELLOW")) %>%
  select(-congruency) %>%
  table

trials2 %>% 
  filter(!(target == "green" & distractor == "PINK")) %>%
  filter(!(target == "yellow" & distractor == "BLACK")) %>%
  select(-congruency) %>%
  table




labels %>% 
  filter(run == "run2") %>%
  filter(!(target == "green" & distractor == "PINK")) %>%
  filter(!(target == "yellow" & distractor == "BLACK")) %>%
  select(target, distractor) %>%
  table


m <- labels %>% 
  filter(run == "run1") %>%
  # filter(!(target == "pink" & distractor == "GREEN")) %>%
  # filter(!(target == "black" & distractor == "YELLOW")) %>%
  select(target, distractor) %>%
  table %>%
  as.matrix


# (rowSums(m) > 1) %*% (t(colSums(m)) > 1) * m


labels <- lapply(labels, as.factor)
labels_a <- data.frame(
  target = paste0(labels$target, "_", labels$congruency),
  distractor = paste0(labels$distractor, "_", labels$congruency),
  stringsAsFactors = TRUE
)

## transform labels to dummy codes:

is_level <- model.matrix(~ . + 0, labels, contrasts.arg = lapply(labels, contrasts, contrasts = FALSE)) == 1
is_level_a <- model.matrix(~ l;. + 0, labels_a, contrasts.arg = lapply(labels_a, contrasts, contrasts = FALSE)) == 1
