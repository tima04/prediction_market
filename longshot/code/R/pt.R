library(data.table)
library(stats4)
library(nleqslv)
library(Rcpp)

source("util.R")
sourceCpp("utils.cpp")

w.prelec <- w_prelec_cpp
    
solve.ps.pt.pw <- function(a, b, t, Rs, epsilon=10e-100) {
    n <- length(Rs)
    fun <- function(xs) { # xs = c(ps, w)
        ps <- xs[1:n]
        w = xs[n+1]
        ## ws is ps^a * (1 - exp(-t*(Rs))) + (1 - ps)^b * (1 - exp(-t))
        rslt <- ps^a * (1 - exp(-t*(Rs)))/t + (1 - ps)^b * (1 - exp(t))/t - w
        return (c(rslt, sum(ps)-1))
    }
    ps <- rep(1/n, n) # initial value
    w <- 0.5
    rslt <- nleqslv(c(ps, w), fun)$x[1:n]
    rslt <- ifelse(rslt <= 0, epsilon, rslt)
    rslt
}

## prelec weighting function
solve.ps.pt.prelec <- function(a, b, t, Rs, epsilon=10e-10) {
    n <- length(Rs)
    fun <- function(xs) { # xs = c(ps, w)
        ps <- xs[1:n]
        w = xs[n+1]
        ## ws is w.prelec(ps, a) * cara(Rs, t) + w.prelec((1 - ps),b) * cara(-1, t) 
        rslt <-  w.prelec(ps, a) * cara(Rs, t) + w.prelec((1 - ps),b) * cara(-1, t) - w
        return (c(rslt, sum(ps)-1))
    }
    ps <- rep(1/n, n) # initial value
    w <- 0.5
    rslt <- nleqslv(c(ps, w), fun)$x[1:n]
    rslt <- ifelse(rslt <= 0, epsilon, rslt)
    rslt
}

## prelec weighting function 5 par
solve.ps.pt.prelec.5par <- function(a1, b1, a2, b2, t, Rs, epsilon=10e-10) {
    n <- length(Rs)
    fun <- function(xs) { # xs = c(ps, w)
        ps <- xs[1:n]
        w = xs[n+1]
        ## ws is w.prelec(ps, a) * cara(Rs, t) + w.prelec((1 - ps),b) * cara(-1, t) 
        rslt <-  w.prelec(ps, a1, b1) * cara(Rs, t) + w.prelec((1 - ps),a2, b2) * cara(-1, t) - w
        return (c(rslt, sum(ps)-1))
    }
    ps <- rep(1/n, n) # initial value
    w <- 0.5
    rslt <- nleqslv(c(ps, w), fun)$x[1:n]
    rslt <- ifelse(rslt <= 0, epsilon, rslt)
    rslt
}


## 6 par (seperate utilities for gain and loss and separate weighting function), u(-1) = lambda
solve.ps.pt.prelec.6par <- function(a1, b1, a2, b2, t, lambda, Rs, epsilon=10e-10) {
    n <- length(Rs)
    fun <- function(xs) { # xs = c(ps, w)
        ps <- xs[1:n]
        w = xs[n+1]
        ## ws is w.prelec(ps, a1, b1) * cara(Rs, t) + w.prelec((1 - ps),a2, b2) * lambda(=u(-1)) 
        rslt <-  w.prelec(ps, a1, b1) * cara(Rs, t) + w.prelec((1 - ps),a2, b2) * lambda - w
        return (c(rslt, sum(ps)-1))
    }
    ps <- rep(1/n, n) # initial value
    w <- 0.5
    rslt <- nleqslv(c(ps, w), fun)$x[1:n]
    rslt <- ifelse(rslt <= 0, epsilon, rslt)
    rslt
}

loglike.pt <- function(a, b, t, dt) { 
    dt[, ps:=solve.ps.pt.pw(a, b, t, winodds), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

## 5 par
loglike.pt.5par <- function(a1,a2,b1,b2,t,dt) { 
    dt[, ps:=solve.ps.pt.prelec.5par(a1,a2,b1,b2,t, winodds), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

## 6 par
loglike.pt.6par <- function(a1,a2,b1,b2,t,lambda,dt) { 
    dt[, ps:=solve.ps.pt.prelec.6par(a1,a2,b1,b2,t,lambda,winodds), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

loglike.pt.6par.parallel <- function(a1,a2,b1,b2,t,lambda,dt) { 

    fun <- function(index) {
        dt <- dt[index]
        dt[, ps:=solve.ps.pt.prelec.6par(a1,a2,b1,b2,t,lambda,winodds), by=raceid]
        rslt <- dt[winner==1, sum(log(ps))]
        rslt
    }
    
    ##dt[ ,utils:=cara_cpp(winodds+1, t)]
    sums <- parLapply(cl, splitIndices(nrow(dt), NNODE), fun=fun)

    rslt <- 0
    for (x in sums)
        rslt <- rslt + x
    rslt
}


dt <- as.data.table(as.data.frame(df.salanie))
##dt <- dt[year==95]
t0 <- date()
out <- mle(minuslogl=function(a,b,t) -loglike.pt(a, b, t, dt), start=list(a=1, b=0.318, t=-0.072),
           method="Nelder-Mead")
## 5 par
## out <- mle(minuslogl=function(a1,b1,a2,b2,t) -loglike.pt(a1,b1,a2,b2,t, dt),
##            start=list(a1=0.5,b1=0.5,a2=0.5,b1=0.5,b2=0.5,t=0.5), method="SANN")
##out <- mle(minuslogl=function(a,b,t) -loglike.pt(a, b, t, dt), start=list(a=0.5, b=0.5, t=0.5) ,method="L-BFGS-B")
t1 = date()
rand <- rbinom(1, 1000, 0.5)
##typ <- "out-pt-prelec-sann"
##typ <- "out-pt-pw95-"
typ <- "out-pt-pw-"
fl <- paste(typ, rand, ".dat", sep="")
out <- list("t0" = t0, "t1" = t1, "out" = out)
save(out, file=fl)

##================================ playground

## out-all-like-507.dat	out-full-like-reg-505.dat   out-pt2-power-480.txt	 out-reg2-514.txt
## out-full2-like-524.dat	out-full-like-reg-505.txt   out-pt2-prelec-531.dat	 out-sal2-510.dat
## out-full-like-465.dat	out-full-like-sal-506.dat   out-pt2-prelec-531.txt	 out-sal2-510.txt
## out-full-like-465.txt	out-full-like-sal-506.txt   out-pt2-prelec-lbfg-502.dat
## out-full-like-476.dat	out-nonmaiden-like-523.dat  out-pt2-prelec-lbfg-502.txt


## load("logs/out-full-like-reg-505.dat")
## out.reg <- out
## load("logs/out-full-like-sal-506.dat")
## out.sal<- out
## load("out-full2-like-524.dat")
## out.pt<- out
## load("logs/out-pt2-power-480.dat")
## out.pt.pow <- out
## load("logs/out-pt2-prelec-531.dat")
## out.pt.pr <- out
## load("logs/out-pt2-prelec-lbfg-502.dat")
## out.pt.pr.lbfg <- out

## summary(out.reg)
## summary(out.sal)
## summary(out.pt)
## summary(out.pt.pow)
## summary(out.pt.pr)
## summary(out.pt.pr.lbfg)
