library(here)
library(data.table)

source(here::here("src", "stroop-rsa-pc.R"))

ttypes <- sort(ttypes)

target <- gsub("[A-Z]", "", ttypes)
distractor <- gsub("[a-z]", "", ttypes)
incongruency <- target != tolower(distractor)
X_target <- tcrossprod(indicator_matrix(target))
X_distractor <- tcrossprod(indicator_matrix(distractor))
X_incongruency <- tcrossprod(indicator_matrix(as.character(incongruency))[, "TRUE"])

colnames(X_target) <- ttypes
colnames(X_distractor) <- ttypes
colnames(X_incongruency) <- ttypes

X_target <- 1 - X_target
X_distractor <- 1 - X_distractor
X_incongruency <- 1 - X_incongruency

image(X_target)
image(X_distractor)
image(X_incongruency)

fwrite(as.data.table(X_target), here("in", "model_target.csv"))
fwrite(as.data.table(X_distractor), here("in", "model_distractor.csv"))
fwrite(as.data.table(X_incongruency), here("in", "model_incongruency.csv"))
