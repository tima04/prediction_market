library(data.table)
library(stats4)
library(nleqslv)
source("util.R")

cara <- function(xs, theta) {
   (1 - exp(-theta*xs))/theta
    ##(1 - exp(-theta*xs))
}

solve.ps <- function(t, gamma, Rs, epsilon=10e-10) {
    n <- length(Rs)
    ts <- t + gamma * seq(0, 1, length=n)
    us <- cara(Rs + 1, ts)
    ## k <- 2 
    ## us[k] <- cara(Rs[k]+1, t + gamma)  
    u.min.1 <- cara(-1, t)
    w <- 1/sum(1/(us - u.min.1)) + u.min.1 
    ps <- (w - u.min.1)/(us - u.min.1)
    return(ps)
}

loglike.eu <- function(t, gamma, dt) { 
    dt[, ps:=solve.ps(t, gamma, winodds), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

dt <- as.data.table(as.data.frame(df))
#out <- mle(minuslogl=function(t) -loglike.eu(t, dt), start=list(t=0.5) ,method="Brent", lower=-1, upper=1)
system.time(out <- mle(minuslogl=function(t, gamma) -loglike.eu(t, gamma, dt),
           start=list(t=0.5, gamma=0) ,method="L-BFGS-B"))


rand <- rbinom(1, 1000, 0.5)
typ <- "out-eu-filter-"
fl <- paste("logs/", typ, rand, ".dat", sep="")
fl.txt <- paste("logs/", typ, rand, ".txt", sep="")
save(out, file=fl)
dput(summary(out), file=fl.txt)

### =================================================================== playground
