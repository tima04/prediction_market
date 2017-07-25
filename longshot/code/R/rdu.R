library(data.table)
library(stats4)
library(nleqslv)
library(Rcpp)
library(parallel)

source("util.R")
sourceCpp("utils.cpp") 

## make clusters, cl is closed at the end of the file
NNODE = detectCores() 
cl <- makeForkCluster(NNODE)
##try(stopCluster(cl), silent = TRUE)

### The utilty function I used is u(x, theta) = cara(x+1, theta), which is in cara family. 

w.prelec <- function(ps, alpha, beta=1) {
    exp(-beta*(-log(ps))^alpha)
}

## dummy.arg not used, just to keep same interface with w.prelec 
w.pow <- function(ps, alpha, dummy.arg=1) {
    ps^alpha
}

w.pow.inv <- function(ps, alpha, dummy.arg=1) {
   ps^(1/alpha) 
}

w.pow.gradient <- function(ps, alpha, dummy.arg=1) {
    alpha * ps^(alpha-1)
}

w.pow.inv.gradient <- function(ps, alpha, dummy.arg=1) {
    w.pow.gradient(ps, 1/alpha)
}

w.prelec.gradient <- function(ps, alpha, beta=1) {
    ((alpha * beta) / ps) * w_prelec_cpp(ps, alpha) * (-log(ps))^(alpha-1)
}

w.prelec.inv.gradient <- function(ps, alpha, beta=1) {
    w.prelec.gradient(ps, 1/alpha, (1/beta)^(1/alpha))
}

## somehow is giving (sligthly) wrong numbers (not matching with EU)
solve.ps.rdu.old <- function(a, utils, prelec = TRUE, epsilon=10e-10) {
    if (prelec)
        wght.fun <- w_prelec_cpp
    else
        wght.fun <- w.pow
    fun <- function(xs) { # xs = c(ps, w)
        ps <- xs[1:n]
        w = xs[n+1]
        ## ws is w.prelec(ps, a) * cara(Rs+1, t)
        wghts <- wght.fun(ps, a)
        rslt <-  wghts * utils - w 
        return (c(rslt, sum(ps)-1))
    }
    n <- length(utils)
    ps <- rep(1/n, n) # initial value
    w <- 0.5
    rslt <- tryCatch(nleqslv(c(ps, w), fun)$x[1:n],
                     error = function(e) {
                         rep(1/n, n)
                         })
    rslt <- ifelse(rslt <= 0, epsilon, rslt)
    rslt
}

## since both prelec and power wghting functions are invertible, inverting them
## simplify solving nonlinear equations, only one equation need to be solved.
solve.ps.rdu <- function(a, beta=1, utils, prelec = TRUE, epsilon=10e-10) {
    if (prelec) {
        w.inv <- w_prelec_inv_cpp
        w.inv.grad <- w.prelec.inv.gradient 
        }
    else {
        w.inv <- w.pow.inv
        w.inv.grad <- w.pow.inv.gradient
        }

    fun <- function(w) { 
        sum(w.inv(w/utils, a, beta)) - 1
    }

    fun.grad <- function(w) {
        sum(w.inv.grad(w/utils, a, beta))
    }

    #w <- tryCatch(nleqslv(1, fun)$x[1],
    w <- tryCatch(nleqslv(1, fun, jac=fun.grad)$x[1],
                     error = function(e) {
                        NA 
                     })

    n <- length(utils)
    if (is.na(w))
        return (rep(1/n, n))
    ps <- w.inv(w/utils, a, beta)
    ps
}


loglike.rdu <- function(a, beta=1, t, dt, prelec=TRUE) { 
    dt[ ,utils:=cara_cpp(winodds+1, t)]
    dt[, ps:=solve.ps.rdu(a, beta=beta, utils=utils, prelec=prelec), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}


## parallel version of loglike 
loglike.rdu.parallel <- function(a, beta, t, dt, prelec=TRUE) { 

    fun <- function(index) {
        dt <- dt[index]
        dt[,ps:=solve.ps.rdu(a, beta, utils, prelec), by=raceid]
        rslt <- dt[winner==1, sum(log(ps))]
        rslt
    }
    
    dt[ ,utils:=cara_cpp(winodds+1, t)]
    sums <- parLapply(cl, splitIndices(nrow(dt), NNODE), fun=fun)

    rslt <- 0
    for (x in sums)
        rslt <- rslt + x
    rslt
}

##dt <- as.data.table(as.data.frame(df))
##dt <- dt[raceid in 1:1000] 
##system.time(

main <- function() {
    dt <- as.data.table(as.data.frame(df.salanie))
    setkey(dt, "raceid")
    t0 <- date()
    ## out <- mle(minuslogl=function(a, beta, t) -loglike.rdu(a, beta, t, dt),
    ##             start=list(a=0.5, beta=0.8, t=-0.02),
    ## ## ##           method = "SANN")
    ##             method = "Nelder-Mead")
    ## method="L-BFGS-B", lower=c(0.4, 0.4, -2), upper=c(2, 2, 2))

    outs <- list()
    for (yr in dt[,unique(year)]) {
        print(yr)
        dt.sub <- dt[year != yr, ]
        outs[yr] <- mle(minuslogl=function(a, beta, t) -loglike.rdu(a, beta, t, dt.sub),
                        start=list(a=0.5, beta=0.8, t=-0.02),
                        method = "Nelder-Mead")
        print(outs[yr])
    }
    t1 = date()
    ## foo <- mcparallel(
    ##     out <- mle(minuslogl=function(a,t) -loglike.rdu2(a, t, dt),
    ##                start=list(a=0.5, t=-0.02) , method="L-BFGS-B",
    ##                lower=c(0.001, -5), upper=c(5, 5))
    ## )
    ## bar <- mccollect(foo)[[1]]
    rand <- rbinom(1, 1000, 0.5)
    typ <- "outs-4-crossvalidation-rdu-"
    ##typ <- "out-rdu-pow-lbfg"
    ##fl <- paste("logs/", typ, rand, ".dat", sep="")
    fl <- paste(typ, rand, ".dat", sep="")
    ##fl.txt <- paste("logs/", typ, rand, ".txt", sep="")
    out = list("t0" = t0, "t1" = t1, "out" = outs)
    save(out, file=fl)
}

try(stopCluster(cl), silent = TRUE)

##===================================== playground

## foo <- function(x) {
##     print(environment())
## }

## dt.sub <- dt[year==86]

## system.time(
## outs <- mle(minuslogl=function(a, beta, t) -loglike.rdu(a, beta, t, dt.sub),
##                 start=list(a=0.5, beta=0.8, t=-0.02),
##             method = "Nelder-Mead")
## )

## system.time(
## outs.p <- mle(minuslogl=function(a, beta, t) -loglike.rdu.parallel(a, beta, t, dt.sub),
##                 start=list(a=0.5, beta=0.8, t=-0.02),
##             method = "Nelder-Mead")
## )
