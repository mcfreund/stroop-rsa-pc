## libraries

library(here)
library(knitr)
library(ggplot2)
library(mfutils)
library(ggseg)
library(ggsegExtra)
library(ggsegSchaefer)
library(ggsegGlasser)
library(broom)
library(parallel)
library(dplyr)
library(tidyr)
library(data.table)

source(here("src", "stroop-rsa-pc.R"))

## settings
theme_set(theme_bw(base_size = 14))

theme_surface <- list(
  theme(
    axis.text = element_blank(), panel.grid = element_blank(), panel.border = element_blank(),
    axis.ticks = element_blank(), legend.position = c(0.5, 0.5), legend.title = element_text(size = 5), 
    legend.background = element_blank(), legend.text = element_text(size = 5), legend.direction = "horizontal",
    legend.key.height = unit(1/5, "cm"), legend.key.width = unit(1/4, "cm"),
    strip.background = element_blank()
  )
)
theme_basic <- list(
  theme(
    axis.text = element_text(size = 6),
    panel.grid = element_blank(), 
    panel.border = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank()
  )
)

colors_model <- c(
  target = "#1b9e77",
  congruency = "#d95f02",
  distractor = "#7570b3"
)

colors_session <- c(
  baseline = "#1f78b4",
  proactive = "#b2df8a",
  reactive = "#a6cee3"
)


## constants

wave <- "wave1"
atlas_roi_nm <- "glasser2016__superparcel"
suffices <- c(
  default = "__nresamp100__lsall_1rpm__divnormed__demeaned",
  centered = "__nresamp100__lsall_1rpm__divnormed__demeaned__centered"
  )

stimsets <- c("bias", "pc50", "biaspc50")
models <- c("target", "distractor", "congruency")
path_data <- here("out", "decoding")
subjdirs <- list.dirs(path_data)[-1]
path_data <- here("out", "decoding")
fnames <- paste0("trainbaseline__", wave, "__svm__", atlas_roi_nm, suffices, "_summarized.RDS")
does_exist <- file.exists(file.path(subjdirs, fnames))
subjdirs <- subjdirs[does_exist]
subjs <- gsub(".*decoding/", "", subjdirs)

subset_nms <- c("bias", "pc50", "biaspc50")
names(ttypes)[names(ttypes) == "all"] <- "biaspc50"

## data
d <- enlist(names(suffices))
for (i in seq_along(fnames)) {
  l <- mclapply(subjdirs, function(x) readRDS(file.path(x, fnames[i])), mc.cores = 20)
  d[[i]] <- rbindlist(l)
}
d <- rbindlist(d, idcol = "fname")
d <- d %>% dcast(... ~ fname, value.var = "evidence")

d[, run := ifelse(trial > max(trial)/2, "run2", "run1"), by = "session"]

cors <- d %>%
  group_by(subj, session, variable, train_stimset, roi, wave) %>%
  summarize(r = (centered/sd(centered) - default/sd(default))) %>%
  group_by(session, variable, train_stimset, roi, wave) %>%
  summarize(r = mean(r))

cors %>%
  ggplot(aes(variable, r, fill = train_stimset)) +
  geom_boxplot() +
  facet_grid(cols = vars(wave, session))


d %>%

  filter(stimulus %in% ttypes$pc50, train_stimset == "pc50") %>%

  group_by(subj, session, variable, roi, wave) %>%
  summarize(b = mean(centered)) %>%
  group_by(session, variable, roi, wave) %>%
  summarize(stat = t.test(b)$statistic) %>%

  ggplot(aes(roi, stat, color = session)) +
  geom_point() +
  geom_line(aes(group = session)) +
  facet_grid(rows = vars(variable))


d %>%

  filter(stimulus %in% ttypes$pc50, train_stimset == "pc50") %>%

  group_by(subj, session, variable, roi, wave) %>%
  summarize(b = mean(default)) %>%
  pivot_wider(names_from = "session", values_from = "b") %>%
  mutate(b = proactive - reactive) %>%
  group_by(variable, roi, wave) %>%
  summarize(stat = t.test(b)$statistic) %>%

  ggplot(aes(roi, stat, color = variable)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  geom_line(aes(group = variable)) +
  facet_grid(rows = vars(variable))
