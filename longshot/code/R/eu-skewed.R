library(data.table)
library(stats4)
library(nleqslv)
source("util.R")

sourceCpp("utils.cpp") 

## utility(-1) = 0
utility.skew <- function(xs, pars) {
    b1 <- pars[1]
    b2 <- pars[2]
    b3 <- pars[3]
    ys <- (b1 - b2 + b3) + b1*xs + b2*(xs^2) + b3*(xs^3)
    ys
} 

solve.ps.eu.skew <- function(utils) {
    w <- 1/sum(1/(utils)) 
    ps <- w / utils 
    return(ps)
}

loglike.eu.skew <- function(pars, dt) { 
    dt[ ,utils:=utility.skew(winodds+1, pars)]
    dt[, ps:=solve.ps.eu.skew(utils), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

dt <- as.data.table(as.data.frame(df.salanie))
system.time(
    out.skew <- mle(
        minuslogl=function(b1, b2, b3) {pars <- c(b1, b2, b3); -loglike.eu.skew(pars, dt)},
        start=list(b1=1, b2=-1, b3=0.5),
        method="Nelder-Mead"))

rand <- rbinom(1, 1000, 0.5)
typ <- "out-eu-filter-"
fl <- paste("logs/", typ, rand, ".dat", sep="")
fl.txt <- paste("logs/", typ, rand, ".txt", sep="")
save(out, file=fl)
dput(summary(out), file=fl.txt)

### =================================================================== playground

