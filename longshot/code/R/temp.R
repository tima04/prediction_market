library(data.table)
library(stats4)
library(nleqslv)
library(Rcpp)
library(parallel)

source("util.R")
sourceCpp("utils.cpp") 

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
solve.ps.rdu <- function(a, utils, prelec = TRUE, epsilon=10e-10) {
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
solve.ps.rdu2 <- function(a, beta=1, utils, prelec = TRUE, epsilon=10e-10) {
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

NNODE = 4 
cl <- makeForkCluster(NNODE)

loglike.rdu <- function(a, beta=1, t, dt, prelec=TRUE) { 
    fun <- function(indices) {
        dt[indices, ps:=solve.ps.rdu2(a, utils, prelec=prelec, beta=beta), by=raceid]
    }
    dt[ ,utils:=cara_cpp(winodds+1, t)]
    dt[, ps:=solve.ps.rdu2(a, utils, prelec=prelec, beta=beta), by=raceid]
    n <- floor(nrow(dt)/2)
    ind1 <- 1:n
    ind2 <- (n+1):nrow(dt) 
    ind <- list(ind1, ind2)
    clusterApply(cl, ind, fun)
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

stopCluster(cl)

##dt <- as.data.table(as.data.frame(df))
##dt <- dt[raceid in 1:1000] 
##system.time(

main <- function() {
    dt <- as.data.table(as.data.frame(df.salanie))
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


##===================================== playground

## load("out-rdu-prelec-neldermead-fast-nonmaiden474.dat")
## out.nm <- out$out
## load("out-rdu-prelec-neldermead-fast-nonmaiden-beta504.dat")
## out.nm.b <- out$out
## load("out-rdu-prelec-neldermead-fast506.dat")
## out.all <- out$out 
## load("out-rdu-prelec-neldermead-fast-salanie512.dat")
## out.sal <- out$out
## load("out-rdu-prelec-neldermead-fast-salanie-beta480.dat")
## out.sal.b <- out$out
## load("out-rdu-prelec-neldermead-fast-maiden-beta498.dat")
## out.m.b <- out$out
## load("out-rdu-prelec-lbfgs-fast-maiden-beta537.dat")
## out.m.b.lb <- out$out
## load("out-rdu-prelec-lbfgs-fast-salanie-beta505.dat")
## out.sal.b.lb <- out$out
## load("out-rdu-prelec-lbfgs-fast-nonmaiden-beta492.dat")
## out.nm.b.lb <- out$out

## ##nm
## a <- 0.1520801
## b <- 0.2770332 

## ## maiden
## a <- 0.2435195 
## b <- 0.1447918 

## a <- 0.8515
## #b <- 1.85
## b <- 1.18569 
## ps <- seq(0,1,0.01)
## plot(w.prelec(ps, a, b)~ps, type="l")
## abline(a=0, b=1)

## ps <- c(0.01,0.1, 0.6)
## w.prelec(ps, a, b)/ps

foo <- function(x) {
    if(x < 0)
        stop("x negative")
    x^2
}
