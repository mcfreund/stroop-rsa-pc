library(colorout)
library(here)
library(data.table)
library(rmarkdown)

source(here("src", "stroop-rsa-pc.R"))

src_dir <- "group_analyses/MCvMI_network_2021-12-16"
ttype_subsets <- c("bias", "pc50")
prewhs <- c("obsbias", "obspc50", "none")
measures <- c("crcor", "cveuc")

params <- expand.grid(ttype_subset = ttype_subsets, prewh = prewhs, measure = measures, stringsAsFactors = FALSE)
for (i in seq_len(nrow(params))) {
    p <- c(params[i, ])
    render_report(name = "MCMI_Schaefer2016Network", src_dir = src_dir, params = p)
    render_report(name = "MIMC_Schaefer2016Network", src_dir = src_dir, params = p)
}
