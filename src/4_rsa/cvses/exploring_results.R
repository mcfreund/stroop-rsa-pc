library(here)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(mfutils)

decoder <- "svm"
atlas_name <- "schaefer2018_17_200"
roi_col <- "parcel"
n_resamples <- 100
glmname <- "lsall_1rpm"
fname <- paste0(decoder, "__", atlas_name, "__", roi_col, "__nresamp", n_resamples, "__", glmname, ".txt")
outd <-fread(here("out", "res_decoding", fname))

outd[, is_farout := farout(distn - mean(distn)), by = c("subj", "roi", "session", "variable", "stimset")]
outd <- outd[is_farout == FALSE]
outd[, distn_s := distn / sd(distn), by = c("subj", "roi", "session", "variable", "stimset")]

d_sum <- outd %>%
    group_by(subj, roi, session, variable, stimset) %>%
    summarize(d = mean(distn_s))

d_stat <- d_sum %>%
    group_by(roi, session, variable, stimset) %>%
    summarize(stat = t.test(d)$statistic, p = wilcox.test(d, alternative = "greater")$p.value) %>%
    group_by(session, variable, stimset) %>%
    mutate(p_fdr = p.adjust(p, "fdr"))

roi_congruency <- intersect(
    d_stat %>% filter(p_fdr < 0.05, session == "proactive", variable == "congruency") %>% pull(roi),
    d_stat %>% filter(p_fdr < 0.05, session == "reactive", variable == "congruency") %>% pull(roi)
)

roi_distractor <- intersect(
    d_stat %>% filter(p_fdr < 0.05, session == "proactive", stimset == "pc50", variable == "distractor") %>% pull(roi),
    d_stat %>% filter(p_fdr < 0.05, session == "reactive", stimset == "pc50", variable == "distractor") %>% pull(roi)
)

roi_target <- intersect(
    d_stat %>% filter(p_fdr < 0.05, session == "proactive", stimset == "pc50", variable == "target") %>% pull(roi),
    d_stat %>% filter(p_fdr < 0.05, session == "reactive", stimset == "pc50", variable == "target") %>% pull(roi)
)
unique(d_stat$roi)

d_sum %>%
    #filter(variable == "target", stimset == "pc50") %>%
    ggplot(aes(roi, d)) +
    stat_summary(fun = mean, geom = "col") +
    facet_grid(rows = vars(variable, stimset), cols = vars(session)) +
    coord_flip()

outd %>%
    filter(stimset != "bias") %>%
    group_by(roi, session, variable, stimset) %>%
    summarize(distn = mean(distn)) %>%
    as.data.table %>%
    dcast(roi + variable + stimset ~ session, value.var = "distn", fun.aggregate = mean) %>%
    ggplot(aes(proactive, reactive)) +
    geom_abline() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point() +
    facet_grid(cols = vars(variable, stimset))


## contrast pro vs rea


outd_contr <- outd %>%
    filter(stimset != "bias") %>%
    dcast(subj + roi + variable + stimset ~ session, value.var = "distn", fun.aggregate = mean) %>%
    mutate(pvr = proactive - reactive) %>%
    group_by(roi, variable, stimset) %>%
    summarize(stat = t.test(pvr)$stat, p = wilcox.test(pvr)$p.value, .groups = "drop")

rois_pfc <- unique(outd_contr$roi)[grep("PFC", unique(outd_contr$roi))]

outd_contr %>%
    ggplot(aes(roi, stat)) +
    geom_col() +
    facet_grid(cols = vars(variable, stimset)) +
    coord_flip()

outd_contr %>%
    filter(variable == "target", roi %in% rois_pfc) %>%
    ggplot(aes(roi, stat)) +
    geom_col() +
    facet_grid(cols = vars(variable, stimset)) +
    coord_flip()

outd_contr %>%
    filter(variable == "distractor", roi %in% rois_pfc) %>%
    ggplot(aes(roi, stat)) +
    geom_col() +
    facet_grid(cols = vars(variable, stimset)) +
    coord_flip()

outd_contr %>%
    filter(variable == "congruency", roi %in% roi_congruency) %>%
    ggplot(aes(roi, stat)) +
    geom_col() +
    facet_grid(cols = vars(variable, stimset)) +
    coord_flip()




outd %>%
    filter(stimset != "bias") %>%
    dcast(subj + roi + variable + stimset ~ session, value.var = "distn", fun.aggregate = mean) %>%
    mutate(pvr = proactive - reactive) %>%
    group_by(roi, variable, stimset) %>%
    summarize(stat = t.test(pvr)$stat, .groups = "drop")


outd_contr %>%
    filter(variable == "target", stimset == "pc50", roi == "R_FEF") %>%
    ggplot(aes(distn, fill = session)) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.5)

outd_contr <- outd_contr %>%
    group_by(roi, variable, stimset) %>%
    mutate(p_fdr = p.adjust(p, "fdr"))

outd_contr %>% filter(p_fdr < 0.05) %>% View






# outd %>%
#     filter(variable == "target", stimset == "pc50", roi == "R_FEF") %>%
#     ggplot(aes(distn, fill = session)) +
#     geom_histogram(bins = 50, position = "identity", alpha = 0.5)

# d_sum %>%
#     filter(variable == "congruency", stimset == "biaspc50", roi == "R_8C") %>%
#     ggplot(aes(d, fill = session)) +
#     geom_histogram(bins = 10, position = "identity", alpha = 0.5)

# outd %>%
#     filter(variable == "congruency", stimset == "biaspc50", roi == "R_8C") %>%
#     group_by(session, subj) %>%
#     ggplot(aes(distn_s, fill = session)) +
#     geom_histogram(bins = 40, position = "identity", alpha = 0.5)





d_sum %>%
    filter(variable == "target", stimset == "pc50", roi %in% c("R_p9-46v")) %>% View


outd %>% filter(variable == "target", stimset == "pc50", roi == "R_FEF", subj == "729254") %>% pull(distn_s) %>% plot()

outd_contr %>%
    filter(variable == "target", roi %in% roi_target, p < 0.05)


outd[, max(distn_s), by = subj]

summary(outd$distn_s)
outd %>%
    ggplot(aes(session, distn)) +
    geom_boxplot()

d_sum <- outd %>%
    group_by(subj, roi, session, variable, stimset) %>%
    mutate(d = mean(distn) / sd(distn))
