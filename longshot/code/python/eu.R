library(data.table)
library(stats4)
library(nleqslv)
source("util.R")

cara <- function(xs, theta) {
    (1 - exp(-theta*xs))/theta
    ##(1 - exp(-theta*xs))
}

solve.ps <- function(t, Rs, epsilon=10e-10) {
    n <- length(Rs)
    us <- cara(Rs + 1, t)
    u.min.1 <- cara(-1, t)
    w <- 1/sum(1/(us - u.min.1)) + u.min.1 
    ps <- (w - u.min.1)/(us - u.min.1)
    return(ps)
}

## solve.ps <- function(t, Rs, epsilon=10e-10) { 
solve.ps.old <- function(t, Rs, epsilon=10e-10) {
    n <- length(Rs)
    fun <- function(xs) { # xs = c(ps, w)
        ps <- xs[1:n]
        w = xs[n+1]
        ## ws is w.prelec(ps, a) * cara(Rs, t) + w.prelec((1 - ps),b) * cara(-1, t) 
        rslt <-  ps * cara(Rs, t) + (1 - ps) * cara(-1, t) - w
        return (c(rslt, sum(ps)-1))
    }
    ps <- rep(1/n, n) # initial value
    w <- 0.5
    rslt <- nleqslv(c(ps, w), fun)$x[1:n]
    rslt <- ifelse(rslt <= 0, epsilon, rslt)
    rslt
}


loglike.eu <- function(t, dt) { 
    dt[, ps:=solve.ps(t, winodds), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
 }

dt <- as.data.table(as.data.frame(df))
out <- mle(minuslogl=function(t) -loglike.eu(t, dt), start=list(t=0.5) ,method="Brent", lower=-1,
           upper=1)
#out <- mle(minuslogl=function(t) -loglike.eu(t, dt), start=list(t=0.5) ,method="L-BFGS-B")


rand <- rbinom(1, 1000, 0.5)
typ <- "out-eu-filter-"
fl <- paste("logs/", typ, rand, ".dat", sep="")
fl.txt <- paste("logs/", typ, rand, ".txt", sep="")
save(out, file=fl)
dput(summary(out), file=fl.txt)

### =================================================================== playground

## ev

