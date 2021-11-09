
## converting unscaled x-y distance to scaled x-y distance

f <- function(d2, lx, ly) (d2 - lx^2 - ly^2)/(lx*ly) + 2
nsim <- 1000
res <- matrix(NA, ncol = 2, nrow = nsim)
for (i in 1:nsim) {
    x <- rnorm(100)
    y <- rnorm(100)
    lx <- sqrt(sum(x^2))
    ly <- sqrt(sum(y^2))
    x_s <- x / lx
    y_s <- y / ly
    d2 <- sum((x - y)^2)
    d2_scale <- sum((x_s - y_s)^2)
    d2_scale_form <- f(d2, lx, ly)
        
    res[i, ] <- c(d2_scale, d2_scale_form)

}
plot(res)
all.equal(res[, 1], res[, 2])

d2_scale


## now for the cross-validated form....

