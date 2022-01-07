library(mfutils)
library(here)
library(data.table)
source("src/stroop-rsa-pc.R")
ssqe <- function(x, y) sum((x - y)^2)

n <- 100  ## trials
p <- 10  ## vertices
B1 <- matrix(rnorm(n*p), nrow = p)
B2 <- matrix(rnorm(n*p), nrow = p)
#B2 <- matrix(rexp(n*p), nrow = p)
S1 <- cov(t(B1))
S2 <- cov(t(B2))
W1 <- solve(S1)
W2 <- solve(S2)
S_bar <- (S1 + S2) / 2
S_bar_inv <- solve(S_bar)
W_bar <- (W1 + W2) / 2

ssqe(W_bar, S_bar_inv)
cor(squareform(W_bar), squareform(S_bar_inv))
plot(squareform(W_bar), squareform(S_bar_inv))



D1 <- t(B1) %*% W1 %*% B2
D2 <- t(B1) %*% W2 %*% B2
D_bar <- (D1 + D2) / 2
D_wbar <- t(B1) %*% W_bar %*% B2
D_sbarinv <- t(B1) %*% S_bar_inv %*% B2

ssqe(D_bar, D_wbar)
ssqe(D_bar, D_sbarinv)
cor(squareform(D_bar), squareform(D_sbarinv))
plot(squareform(D_bar), squareform(D_sbarinv))

