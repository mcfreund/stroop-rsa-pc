
# library(colorout)
# library(here)
# library(dplyr)
# library(tidyr)
# library(data.table)
# library(gifti)
# library(abind)



library(colorout)
library(here)
source(here("src", "stroop-rsa-pc.R"))

atlas <- read_atlas()
files <- c(
    "/data/nil-external/ccp/freund/stroop-rsa-pc/out/glms/115825/wave1/RESULTS/Stroop/reactive_lsall_1rpm/STATS_115825_run1_L_REML.func.gii",
    "/data/nil-external/ccp/freund/stroop-rsa-pc/out/glms/115825/wave1/RESULTS/Stroop/reactive_lsall_1rpm/STATS_115825_run1_R_REML.func.gii",
    "/data/nil-external/ccp/freund/stroop-rsa-pc/out/glms/115825/wave1/RESULTS/Stroop/reactive_lsall_1rpm/STATS_115825_run2_L_REML.func.gii",
    "/data/nil-external/ccp/freund/stroop-rsa-pc/out/glms/115825/wave1/RESULTS/Stroop/reactive_lsall_1rpm/STATS_115825_run2_R_REML.func.gii"
)
l <- lapply(files, gifti::read_gifti)
gii <- l[[1]]

new_parcellated_data()
str(extract_labels(gii))
str(extract_data(gii, pattern = "_Coef"))
str(concat_hemis(l[1:2], pattern = "_Coef"))
str(parcellate_data(concat_hemis(l[1:2], pattern = "_Coef"), atlas))
str(parcellate_data(concat_hemis(l[1:2], pattern = "_Coef"), atlas))
x <- giftis_to_parcellated_image(
    giftis = l, folds = c("run1", "run1", "run2", "run2"), hemis = c("L", "R", "L" , "R"), atlas = atlas,
    glm_name = "lsall_1rpm", roi_set = "Schaefer2018_control", prewhitened = "none", 
    shrinkage_var = numeric(),
    shrinkage_cov = numeric(),
    subject = "",
    wave = "",
    session = "",
    task = "",
    pattern = "_Coef"
)
x
lapply(x$data, class)
x$bad_vertices


file.exists(construct_filename_gifti("132017", "wave1", "baseline", "run1", "lsall_1rpm", "L"))
file.exists(construct_filename_gifti("132017", "wave1", c("baseline", "baseline"), "run1", "lsall_1rpm", "L"))

file.exists(construct_filename_h5("pimage-B", "lsall_1rpm", "Schaefer2018_control", "none"))
file.exists(construct_filename_h5("pimage-B", "lsall_1rpm", "Schaefer2018_control", c("none", "none")))

create_nested_group()


a <- matrix(rnorm(1E6), ncol = 200)
nms <- write_dset(
    a, dset_prefix = "coefs",
    subject = "132017", wave = "wave1", session = "baseline", run = "run1", 
    roiset = "Schaefer2018_control", roi = "Vis",
    glmname = "lsall-1rpm", prewh = "none"
    )
write_links(list(nms))
h5ls(here::here("out", "parcellated", "master.h5"))
