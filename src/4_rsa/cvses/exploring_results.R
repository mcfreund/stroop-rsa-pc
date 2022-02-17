library(here)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(mfutils)

theme_set(theme_bw())

decoder <- "svm"
atlas_name <- "glasser2016"
roi_col <- "parcel"
n_resamples <- 100
glmname <- "lsall_1rpm"
fname <- paste0(decoder, "__", atlas_name, "__", roi_col, "__nresamp", n_resamples, "__", glmname, "__scaled.txt")
outd <- fread(here("out", "res_decoding", fname))

outd[, is_farout := farout(distn - mean(distn)), by = c("subj", "roi", "session", "variable", "stimset")]
outd <- outd[is_farout == FALSE]
outd[, distn_s := distn / sd(distn), by = c("subj", "roi", "session", "variable", "stimset")]

d_sum <- outd %>%
    group_by(subj, roi, session, variable, stimset) %>%
    summarize(d = mean(distn_s))

outd_contr <- outd %>%
    filter(stimset != "bias") %>%
    as.data.table %>%
    dcast(subj + roi + variable + stimset ~ session, value.var = "distn_s", fun.aggregate = mean) %>%
    mutate(pr_sum = proactive + reactive, pr_diff = proactive - reactive)

d_stat <- outd_contr %>%
    group_by(roi, variable) %>%
    summarize(stat = t.test(pr_sum)$statistic, p = wilcox.test(pr_sum, alternative = "greater")$p.value) %>%
    group_by(variable) %>%
    mutate(p_fdr = p.adjust(p, "fdr"))

roi_congruency <- d_stat %>% filter(p_fdr < 0.05, variable == "congruency") %>% pull(roi)
roi_distractor <- d_stat %>% filter(p_fdr < 0.05, variable == "distractor") %>% pull(roi)
roi_target <- d_stat %>% filter(p_fdr < 0.05, variable == "target") %>% pull(roi)







d_stat <- d_sum %>%
    group_by(roi, session, variable, stimset) %>%
    summarize(stat = t.test(d)$statistic, p = wilcox.test(d, alternative = "greater")$p.value) %>%
    group_by(session, variable, stimset) %>%
    mutate(p_fdr = p.adjust(p, "fdr"))

roi_target <- union(
    d_stat %>% filter(p_fdr < 0.05, session == "proactive", stimset == "pc50", variable == "target") %>% pull(roi),
    d_stat %>% filter(p_fdr < 0.05, session == "reactive", stimset == "pc50", variable == "target") %>% pull(roi)
)
d_stat


# outd_contr <- outd %>%
#     filter(stimset != "bias") %>%
#     as.data.table %>%
#     dcast(subj + roi + variable + stimset ~ session, value.var = "distn_s", fun.aggregate = mean) %>%
#     mutate(pr_sum = proactive + reactive, pr_diff = proactive - reactive)

# d_stat <- outd_contr %>%
#     group_by(roi, variable) %>%
#     summarize(stat = t.test(pr_sum)$statistic, p = wilcox.test(pr_sum, alternative = "greater")$p.value) %>%
#     group_by(variable) %>%
#     mutate(p_fdr = p.adjust(p, "fdr"))
# roi_congruency <- d_stat %>% filter(p_fdr < 0.05, variable == "congruency") %>% pull(roi)
# roi_distractor <- d_stat %>% filter(p_fdr < 0.05, variable == "distractor") %>% pull(roi)
# roi_target <- d_stat %>% filter(p_fdr < 0.05, variable == "target") %>% pull(roi)

# outd_contr %>%
#     filter(roi %in% roi_congruency) %>%
#     ggplot(aes(roi, pr_diff)) +
#     stat_summary(fun = mean, geom = "col") +
#     stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0, color = "grey30") +
#     coord_flip()


outd_contr %>%
    filter(roi %in% roi_congruency, variable == "congruency") %>%
    ggplot(aes(roi, pr_diff)) +
    stat_summary(fun = mean, geom = "col") +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0, color = "grey30") +
    coord_flip()


outd_contr %>%
    filter(roi %in% roi_target, variable == "target") %>%
    group_by(roi, variable) %>%
    summarize(stat = t.test(pr_diff)$statistic, p = wilcox.test(pr_diff)$p.value) %>%
    group_by(variable) %>%
    mutate(p_fdr = p.adjust(p, "fdr")) %>%
    filter(p < 0.05)



#d_stat <- 
# d_sum %>%
#     filter(stimset != "bias") %>%
#     pivot_wider(id_cols = c("subj", "roi", "variable"), names_from = "session", values_from = "d") %>%
#     mutate(ppv)
#     group_by(roi, session, variable, stimset) %>%
#     summarize(stat = t.test(d)$statistic, p = wilcox.test(d, alternative = "greater")$p.value) %>%
#     group_by(session, variable, stimset) %>%
#     mutate(p_fdr = p.adjust(p, "holm"))

d_stat <- d_sum %>%
    group_by(roi, session, variable, stimset) %>%
    summarize(stat = t.test(d)$statistic, p = wilcox.test(d, alternative = "greater")$p.value) %>%
    group_by(session, variable, stimset) %>%
    mutate(p_fdr = p.adjust(p, "fdr"))

roi_congruency <- intersect(
    d_stat %>% filter(p_fdr < 0.05, session == "proactive", variable == "congruency") %>% pull(roi),
    d_stat %>% filter(p_fdr < 0.05, session == "reactive", variable == "congruency") %>% pull(roi)
)

# outd_contr <- d_stat %>%
#     filter(stimset != "bias") %>%
#     as.data.table %>%
#     dcast(subj + roi + variable + stimset ~ session, value.var = "distn_s", fun.aggregate = mean) %>%
#     mutate(pr_sum = proactive + reactive, pr_diff = proactive - reactive)

d_stat_contr <- outd_contr %>%
    group_by(roi, variable) %>%
    summarize(stat = t.test(pr_diff)$statistic, p = wilcox.test(pr_diff, alternative = "greater")$p.value) %>%
    group_by(variable) %>%
    mutate(p_fdr = p.adjust(p, "fdr"))

d_stat_contr %>%
    filter(roi %in% roi_congruency, variable == "target") %>% View


# roi_distractor <- union(
#     d_stat %>% filter(p_fdr < 0.05, session == "proactive", stimset == "pc50", variable == "distractor") %>% pull(roi),
#     d_stat %>% filter(p_fdr < 0.05, session == "reactive", stimset == "pc50", variable == "distractor") %>% pull(roi)
# )

# roi_target <- union(
#     d_stat %>% filter(p_fdr < 0.05, session == "proactive", stimset == "pc50", variable == "target") %>% pull(roi),
#     d_stat %>% filter(p_fdr < 0.05, session == "reactive", stimset == "pc50", variable == "target") %>% pull(roi)
# )
unique(d_stat$roi)

md <- list(
  core       = combo_paste(c("L_", "R_"), c("p9-46v", "a9-46v", "i6-8", "AVI", "8C", "IFJp", "IP2", "IP1", "PFm", "8BM", "SCEF")),
  extended   = combo_paste(c("L_", "R_"), c(
    "a9-46v", "p10p", "a10p", "11l", "a47r", "p47r", "FOP5", "AVI", "p9-46v", "8C", "IFJp", "6r", "s6-8", "i6-8",
    "SCEF", "8BM", "a32pr", "d32",
    "TE1m", "TE1p",
    "AIP", "IP2", "LIPd", "MIP", "IP1", "PGs", "PFm", "POS2"
  ))
)


d_sum %>%
    filter(stimset != "bias") %>%
    ggplot(aes(roi, d, fill = session)) +
    stat_summary(fun = mean, geom = "col", position = position_dodge()) +
    facet_grid(rows = vars(variable, stimset))

d_sum %>%
    filter(stimset != "bias", roi %in% md$core) %>%
    ggplot(aes(roi, d, fill = session)) +
    stat_summary(fun = mean, geom = "col", position = position_dodge(width = 1)) +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0, position = position_dodge(width = 1), color = "grey30") +
    facet_grid(cols = vars(variable, stimset)) +
    scale_fill_viridis_d() +
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


# outd_contr <- outd %>%
#     filter(stimset != "bias") %>%
#     dcast(subj + roi + variable + stimset ~ session, value.var = "distn", fun.aggregate = mean) %>%
#     mutate(pvr = proactive - reactive) %>%
#     group_by(roi, variable, stimset) %>%
#     summarize(stat = t.test(pvr)$stat, p = wilcox.test(pvr)$p.value, .groups = "drop")

outd_contr <- outd %>%
    filter(stimset != "bias") %>%
    dcast(subj + roi + variable + stimset ~ session, value.var = "distn", fun.aggregate = mean) %>%
    mutate(pvr = proactive - reactive)


#rois_pfc <- unique(outd_contr$roi)[grep("PFC", unique(outd_contr$roi))]

outd_contr %>%
    filter(roi %in% md$core) %>%
    ggplot(aes(roi, pvr)) +
    stat_summary(fun = mean, geom = "col") +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0, color = "grey30") +
    facet_grid(cols = vars(variable)) +
    coord_flip()


outd_contr %>%
    filter(roi %in% md$extended) %>%
    group_by(subj, variable) %>%
    summarize(pvr = mean(pvr)) %>%
    ggplot(aes(variable, pvr)) +
    stat_summary(fun = mean, geom = "col") +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0, color = "grey30")




outd_contr %>%
    filter(roi %in% md$core) %>%
    ggplot(aes(roi, stat)) +
    geom_col() +
        stat_summary(fun = mean, geom = "col", position = position_dodge(width = 1)) +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0, position = position_dodge(width = 1), color = "grey30") +

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


outd %>%
    filter(variable == "target", stimset == "pc50", roi == "Cingulo-Opercular") %>%
    group_by(subj, session) %>% summarize(distn = mean(distn)) %>%
    ggplot(aes(distn, fill = session)) +
    geom_histogram(bins = 15, position = "identity", alpha = 0.5)

outd %>%
    filter(variable == "congruency", stimset == "biaspc50", roi == "Language") %>%
    #group_by(subj, session) %>% summarize(distn = mean(distn_s)) %>%
    ggplot(aes(distn_s, fill = session)) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
    scale_fill_viridis_d()



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
