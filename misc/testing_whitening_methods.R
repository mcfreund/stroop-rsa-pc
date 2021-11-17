n <- 100#1E6  ## for numerical accuracy
p <- 2
X <- matrix(rnorm(n*p, sd = 10), ncol = p)
S <- cov(X)

W2 <- solve(S)
W_solve_pracma <- pracma::sqrtm(W2)$B
W_pracma <- pracma::sqrtm(S)$Binv


S_svd <- svd(S)
D <- diag(1/S_svd$d)
U <- svd(S)$u
W_svd_pracma <- tcrossprod(U %*% D, U)  ## pracma help
W_svd_esl <- tcrossprod(D, U)  ## from elements of stat learning


sum(colSums((W_solve_pracma - W_pracma)^2))
sum(colSums((W_svd_esl - W_svd_pracma)^2))  ## equal within numerical error (?)
sum(colSums((W_solve_pracma - W_svd_pracma)^2))  ## not equal

cov(X)
cov(X %*% W_pracma)
cov(X %*% W_solve_pracma)
cov(X %*% W_svd_pracma)
cov(X %*% W_svd_esl)


## pracma is best method.
## other analysis suggested pracma::sqrtm()$Binv is faster than expm::sqrtm(solve(S)), too.