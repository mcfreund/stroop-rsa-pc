library(here)
library(ggplot2)
library(viridis)
library(colorspace)
library(grid)
library(gridExtra)
library(cowplot)
library(reshape2)
library(magrittr)
library(dplyr)
library(purrr)
library(tidyr)
library(mikeutils)

theme_set(theme_bw(base_size = 14))
# colors.model <- c(orange = "#d95f02", green = "#1b9e77", purple = "#7570b3")
text.size <- 3.5



## create omnibus models ----

colors_bias <- c("blue", "purple", "red", "white")
colors_pc50 <- c("black", "green", "pink", "yellow")
words_bias <- toupper(colors_bias)
words_pc50 <- toupper(colors_pc50)
labels_bias <- expand.grid(list(words_bias, colors_bias))
labels_pc50 <- expand.grid(list(words_pc50, colors_pc50))
names(labels_bias) <- c("word", "color")
names(labels_pc50) <- c("word", "color")
labels <- rbind(labels_bias, labels_pc50)
labels$color <- as.character(labels$color)  ## for plotting
is.inc <- labels$color != tolower(labels$word)
labels$congruency <- ifelse(is.inc, "I", "C")
labels$color_agg <- paste0(labels$color, "_", labels$congruency)
labels$word_agg <- paste0(labels$word, "_", labels$congruency)

labels_f <- lapply(labels, as.factor)
is_level <- model.matrix(~ . + 0, labels_f, contrasts.arg = lapply(labels_f, contrasts, contrasts = FALSE)) == 1

empty <- diag(nrow(labels))
colnames(empty) <- apply(labels, 1, paste, collapse = "_")
rownames(empty) <- apply(labels, 1, paste, collapse = "_")

hue <- empty  ## hue
wor <- empty  ## word
inc <- empty

for (col.i in grep("color", colnames(is_level))) hue[is_level[, col.i], is_level[, col.i]] <- 1
for (col.i in grep("word", colnames(is_level))) wor[is_level[, col.i], is_level[, col.i]] <- 1

diag(inc) <- 0
inc[is.inc, is.inc] <- 1



## now for aggregated models

word_agg <- is_level[, grep("word_agg", colnames(is_level))]
colnames(word_agg) <- gsub("word_agg", "", colnames(word_agg))
word_agg <- sweep(word_agg, 2, colSums(word_agg), "/")

color_agg <- is_level[, grep("color_agg", colnames(is_level))]
colnames(color_agg) <- gsub("color_agg", "", colnames(color_agg))
color_agg <- sweep(color_agg, 2, colSums(color_agg), "/")

inc_a <- t(word_agg) %*% inc %*% word_agg
wor_a <- t(word_agg) %*% wor %*% word_agg
hue_a <- t(color_agg) %*% hue %*% color_agg

hue_wordagg <- t(word_agg) %*% hue %*% word_agg
word_hueagg <- t(color_agg) %*% wor %*% color_agg

cor(
  hue_wordagg[lower.tri(hue_wordagg)],
  wor_a[lower.tri(wor_a)]
)
cor(
  hue_wordagg[lower.tri(hue_wordagg)],
  inc_a[lower.tri(wor_a)]
)



##  plot ----


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
  
  cbind(labels, setNames(as.data.frame(c(R)), val.name))
  
}



plotmat <- function(x, title, axis.text.size = 2, strip.text.size = 1.5) {
  
  cols <- ifelse(as.character(labels$color) == "white", "grey50", as.character(labels$color))
  
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
      scale_x_discrete(labels = labels$word) +
      scale_y_discrete(labels = rev(as.character(labels$word))) +
      
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(
        color = cols, face = "bold", size = rel(axis.text.size),
        angle = 90, hjust = 1, vjust = 0.5
        ),
        axis.text.y = element_text(color = rev(cols), face = "bold", size = rel(axis.text.size)),
        axis.title = element_blank(),
        plot.title = element_text(size = rel(2)),
        panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = rel(strip.text.size))
      ) +
      
      labs(title = title)
  
  
}

p_hue <- hue %>%
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1")) %>%
  separate(v2, c("word2", "color2")) %>%
  plotmat("target", axis.text.size = 1)

p_wor <- wor %>%
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1")) %>%
  separate(v2, c("word2", "color2")) %>%
  plotmat("distractor", axis.text.size = 1)

p_inc <- inc %>%
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1")) %>%
  separate(v2, c("word2", "color2")) %>%
  plotmat("incongruency", axis.text.size = 1)


p <- arrangeGrob(
  p_hue,
  p_wor,
  p_inc,
  nrow = 1
  # layout_matrix = rbind(c(1:3, NA), 4:7)
)

ggsave(
  "C:/Users/mcf/Box/global/docs/papers/diss_proposal/figs/a1_models_omnibus.pdf",
  p,
  width = 4*3, height = 4, unit = "in", dev = "pdf"
)



p_hue_pc50 <- hue %>%
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1")) %>%
  separate(v2, c("word2", "color2")) %>%
  filter(color1 %in% colors_pc50, color2 %in% colors_pc50) %>%
  plotmat("target", axis.text.size = 1)

p_wor_pc50 <- wor %>%
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1")) %>%
  separate(v2, c("word2", "color2")) %>%
  filter(color1 %in% colors_pc50, color2 %in% colors_pc50) %>%
  plotmat("distractor", axis.text.size = 1)

p_inc_pc50 <- inc %>%
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1")) %>%
  separate(v2, c("word2", "color2")) %>%
  filter(color1 %in% colors_pc50, color2 %in% colors_pc50) %>%
  plotmat("incongruency", axis.text.size = 1)


p_pc50 <- arrangeGrob(
  p_hue_pc50,
  p_wor_pc50,
  p_inc_pc50,
  nrow = 1
  # layout_matrix = rbind(c(1:3, NA), 4:7)
)

ggsave(
  "C:/Users/mcf/Box/global/docs/papers/diss_proposal/figs/a1_models_pc50.pdf",
  p_pc50,
  width = 4*3, height = 4, unit = "in", dev = "pdf"
)


## for aggregated models


plotmat_agg <- function(x, title, axis.text.size = 2, strip.text.size = 1.5) {
  
  x %>%
    
    mutate(
      
      x = paste0(word1),
      y = paste0(word2),
      
      x = factor(x, levels = sort(unique(x))),
      y = factor(y, levels = rev(sort(unique(y))))
      
    ) %>%
    
    ggplot(aes(x, y, fill = value)) +
    geom_raster() +
    
    scale_fill_viridis_c(option = "magma") +
    # scale_x_discrete(labels = labels$word) +
    # scale_y_discrete(labels = rev(as.character(labels$word))) +
    
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_text(
        face = "bold", size = rel(axis.text.size),
        angle = 90, hjust = 1, vjust = 0.5
      ),
      axis.text.y = element_text(face = "bold", size = rel(axis.text.size)),
      axis.title = element_blank(),
      plot.title = element_text(size = rel(2)),
      panel.spacing = unit(0, "lines"), 
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = rel(strip.text.size))
    ) +
    
    labs(title = title) +
    
    facet_grid(
      vars(congr1), vars(congr2), 
      switch = "both",
      labeller = labeller(congr1 = c(C = "Congr.", I = "Incon."), congr2 = c(C = "Congr.", I = "Incon."))
    )
  
  
}


p_distractor_agg <- hue_a %>%
  symmat4ggplot %>%
  separate(v1, c("word1", "congr1")) %>%
  separate(v2, c("word2", "congr2")) %>%
  plotmat_agg("target/distractor", axis.text.size = 1)

p_incongruent_agg <- inc_a %>%
  symmat4ggplot %>%
  separate(v1, c("word1", "congr1")) %>%
  separate(v2, c("word2", "congr2")) %>%
  plotmat_agg("incongruency", axis.text.size = 1)


inds <- grep(paste0(colors_pc50, collapse = "|"), colnames(hue_a))
p_distractor_agg_pc50 <- hue_a[inds, inds] %>%
  symmat4ggplot %>%
  separate(v1, c("word1", "congr1")) %>%
  separate(v2, c("word2", "congr2")) %>%
  plotmat_agg("target/distractor", axis.text.size = 1)

p_incongruent_agg_pc50 <- inc_a[inds, inds] %>%
  symmat4ggplot %>%
  separate(v1, c("word1", "congr1")) %>%
  separate(v2, c("word2", "congr2")) %>%
  plotmat_agg("incongruency", axis.text.size = 1)



# p_nuisance_agg <-
#   hue_wordagg %>%
#   symmat4ggplot %>%
#   separate(v1, c("word1", "congr1")) %>%
#   separate(v2, c("word2", "congr2")) %>%
#   plotmat_agg("nuisance", axis.text.size = 1)

p_agg <- arrangeGrob(
  p_distractor_agg,
  p_incongruent_agg,
  p_distractor_agg_pc50,
  p_incongruent_agg_pc50,
  # p_nuisance_agg,
  nrow = 2
  # layout_matrix = rbind(c(1:3, NA), 4:7)
)

ggsave(
  "C:/Users/mcf/Box/global/docs/papers/diss_proposal/figs/a1_models_agg.pdf",
  p_agg,
  width = 4*3, height = 4*3, unit = "in", dev = "pdf"
)


# p_agg_pc50 <- arrangeGrob(
#   p_distractor_agg_pc50,
#   p_incongruent_agg_pc50,
#   # p_nuisance_agg_pc50,
#   nrow = 1
#   # layout_matrix = rbind(c(1:3, NA), 4:7)
# )
# 
# ggsave(
#   "C:/Users/mcf/Box/global/docs/papers/diss_proposal/figs/a1_models_pc50_agg.pdf",
#   p_agg_pc50,
#   width = 4*3, height = 4, unit = "in", dev = "pdf"
# )





## VIFs

vif <- function(x) {
  vifs <- car::vif(lm(1:nrow(x) ~ ., as.data.frame(x)))
  data.frame(vifs, model = names(vifs))
}

X_full <- cbind(
  target = hue[lower.tri(hue)],
  distractor = wor[lower.tri(wor)],
  incongruency = inc[lower.tri(inc)]
)

inds_bias <- 1:16
inds_pc50 <- 17:32



lt <- lower.tri(diag(inds_bias))

X_bias <- cbind(
  target = hue[inds_bias, inds_bias][lt],
  distractor = wor[inds_bias, inds_bias][lt],
  incongruency = inc[inds_bias, inds_bias][lt]
)

X_pc50 <- cbind(
  target = hue[inds_pc50, inds_pc50][lt],
  distractor = wor[inds_pc50, inds_pc50][lt],
  incongruency = inc[inds_pc50, inds_pc50][lt]
)


take_out <- which(grepl("blue|red", colnames(hue)) & !grepl("BLUE_blue|RED_red", colnames(hue)))
inds_bias_rea <- setdiff(inds_bias, take_out)

lt_rea <- lower.tri(diag(inds_bias_rea))

X_bias_rea <- cbind(
  target = hue[inds_bias_rea, inds_bias_rea][lt_rea],
  distractor = wor[inds_bias_rea, inds_bias_rea][lt_rea],
  incongruency = inc[inds_bias_rea, inds_bias_rea][lt_rea]
)


vifs <- list(omnibus = X_full, `pc50/bias` = X_pc50, bias_rea = X_bias_rea) %>% 
  lapply(vif) %>% 
  bind_rows(.id = "stimulus_set")

vifs %>%
  
  arrange(stimulus_set, model) %>%
  
  ggplot(aes(vifs, model, fill = stimulus_set)) +
  geom_col(width = 0.5, position = position_dodge()) +
  
  scale_fill_viridis_d() +
  
  theme(legend.position = c(0.5, 1/3)) +
  labs(
    x = "variance inflation factor", fill = "stimulus set"
  )



cbind(
  inc_a[inds, inds]
)


X_bias_rea <- cbind(
  target = hue[inds_bias_rea, inds_bias_rea][lt_rea],
  distractor = wor[inds_bias_rea, inds_bias_rea][lt_rea],
  incongruency = inc[inds_bias_rea, inds_bias_rea][lt_rea]
)


vifs <- list(omnibus = X_full, `pc50/bias` = X_pc50, bias_rea = X_bias_rea) %>% 
  lapply(vif) %>% 
  bind_rows(.id = "stimulus_set")























p <- arrangeGrob(
  p.rule,
  p.inc,
  p.target,
  p.stimcol,
  p.stimwor,
  p.distr,
  p.response,
  layout_matrix = rbind(c(1:3, NA), 4:7)
)

ggsave(
  "C:/Users/mcf/Box/global/apps/r21_stroop-rsa-aging/figs/drafts/allmods_short.pdf",
  p,
  width = 4*4, height = 4*2, unit = "in", dev = "pdf"
)



ggsave(
  "C:/Users/mcf/Box/global/apps/r21_stroop-rsa-aging/figs/drafts/legend.pdf",
  p.rule + theme(legend.position = "right"),
  width = 4, height = 4, unit = "in", dev = "pdf"
)


## data and RDM ----


n.features <- 64
n.conditions <- nrow(labels)
nms <- apply(labels, 1, paste0, collapse = "_")

b.tar <- 0.2
b.rul <- 0.4

S <- tar*b.tar + rul*b.rul
diag(S) <- 1
d <- MASS::mvrnorm(n.features, mu = rnorm(n.conditions), Sigma = S)
colnames(d) <- nms

p_rdm <- as.matrix(dist(t(scale(d)))) %>%
  symmat4ggplot %>%
  separate(v1, c("word1", "color1", "rule1")) %>%
  separate(v2, c("word2", "color2", "rule2")) %>%
  plotmat("") +
  facet_grid(
    vars(rule1), vars(rule2), 
    switch = "both",
    labeller = labeller(rule1 = rule_labs, rule2 = rule_labs)
  )


ggsave(
  "C:/Users/mcf/Box/global/apps/r21_stroop-rsa-aging/figs/drafts/rdm.pdf",
  p_rdm,
  width = 4, height = 4, unit = "in", dev = "pdf"
)


for (ii in 1:ncol(d)) {
  
  d_ii <- d[, ii, drop = FALSE] %>%
    melt %>%
    ggplot(aes(Var1, Var2, fill = value)) +
    geom_raster() +
    scale_fill_viridis_c(option = "magma") +
    
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      plot.title = element_blank(),
      panel.border = element_blank(),
      panel.spacing = unit(0, "lines"), 
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.ticks = element_blank()
    )
  
  
  
  
  ggsave(
    paste0("C:/Users/mcf/Box/global/apps/r21_stroop-rsa-aging/figs/drafts/d", ii, ".pdf"),
    d_ii,
    width = 4, height = 1/2, unit = "in", dev = "pdf"
  )
  
}





## single trial ----


p.rule.single <- 
  
  (1 - rul["BLUE_blue_CN", , drop = FALSE]) %>%
  
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1", "rule1")) %>%
  separate(v2, c("word2", "color2", "rule2")) %>%
  
  plotmat("rule") +
  
  facet_grid(
    vars(rule2), vars(rule1), 
    switch = "both",
    labeller = labeller(rule1 = rule_labs, rule2 = rule_labs)
  )



p.inc.single <- 
  
  (1 - inc["BLUE_blue_CN", , drop = FALSE]) %>%
  
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1", "rule1")) %>%
  separate(v2, c("word2", "color2", "rule2")) %>%
  
  plotmat("incongruency") +
  
  facet_grid(
    vars(rule2), vars(rule1), 
    switch = "both",
    labeller = labeller(rule1 = rule_labs, rule2 = rule_labs)
  )


p.dis.single <- 
  
  (1 - dis["BLUE_blue_CN", , drop = FALSE]) %>%
  
  symmat4ggplot %>%
  separate(v1, c("word1", "color1", "rule1")) %>%
  separate(v2, c("word2", "color2", "rule2")) %>%
  
  plotmat("distractor") +
  
  facet_grid(
    vars(rule2), vars(rule1), 
    switch = "both",
    labeller = labeller(rule1 = rule_labs, rule2 = rule_labs)
  )



# triali <- inc["GREEN_orange_CN", , drop = FALSE]
# triali[] <- c(rexp(length(c(triali))/2), rexp(length(c(triali))/2)+1)


as.matrix(dist(t(scale(d)))) %>%
  symmat4ggplot %>%
  filter(v1 == "BLUE_blue_CN") %>%
  separate(v1, c("word1", "color1", "rule1")) %>%
  separate(v2, c("word2", "color2", "rule2")) %>%
  plotmat("") +
  facet_grid(
    vars(rule1), vars(rule2), 
    switch = "both",
    labeller = labeller(rule1 = rule_labs, rule2 = rule_labs)
  )



p.triali <- 
  
  as.matrix(dist(t(scale(d)))) %>%
  symmat4ggplot %>%
  filter(v1 == "BLUE_blue_CN") %>%
  separate(v1, c("word1", "color1", "rule1")) %>%
  separate(v2, c("word2", "color2", "rule2")) %>%
  
  plotmat("") +
  
  facet_grid(
    vars(rule2), vars(rule1), 
    switch = "both",
    labeller = labeller(rule1 = rule_labs, rule2 = rule_labs)
  )


p.single <- arrangeGrob(
  p.triali,
  p.rule.single,
  p.inc.single,
  p.dis.single,
  ncol = 4
)


ggsave(
  "C:/Users/mcf/Box/global/apps/r21_stroop-rsa-aging/figs/drafts/singlemod.pdf",
  p.single,
  width = 5.75, height = 4, unit = "in", dev = "pdf"
)

