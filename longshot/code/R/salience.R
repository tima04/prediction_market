library(data.table)
library(stats4)
library(nleqslv)
library(Rcpp)
library(parallel)

source("util.R")
sourceCpp("utils.cpp") 

### The utilty function I used is u(x, theta) = cara(x+1, theta), which is in cara family. 
solve.ps.sal <- function(delta, utils, epsilon=10e-10) {
    n <- length(utils)
    #deltas <- pow_seq(delta, n) ## (delta, delta^2,..,delta^n)
    deltas <- delta^(1:n) ## (delta, delta^2,..,delta^n)
    ps <- deltas/utils
    total <- sum(ps) ## normalizing
    ps/total
}

## sig.hat(x_i, f(xs_minus.i)) = sig(x_i, mean(xs_minus.i))
solve.ps.sal.mean <- function(delta, utils, sal.rank, epsilon=10e-10) {
    n <- length(utils)
    deltas <- pow_seq(delta, n) ## (delta, delta^2,..,delta^n)
    deltas <- delta^sal.rank
    ps <- deltas/utils
    total <- sum(ps) ## normalizing
    ps/total
}

get.sal.rank <- function(rs) {
    n <- length(rs)
    rslt <- c()
    for (i in 1:n) {
        others <- (rs - (n-2))/(n-1)
        others[i] <- rs[i]
        rslt <- c(rslt, match(others[i], sort(others, decreasing = TRUE)))
    }
    rslt
}

loglike.sal <- function(delta, t, dt) { 
    dt[ ,utils:=cara_cpp(winodds+1, t)]
    dt[, ps:=solve.ps.sal(delta, utils), by=raceid]
    ##dt[, ps:=solve.ps.sal.mean(delta, utils, sal.rank), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

## dt <- as.data.table(as.data.frame(df.salanie))
## setkey(dt,"raceid", "winodds")
## system.time(dt[,sal.rank := get.sal.rank(winodds), by=raceid])

## ## dividing racing into two groups, with heavy-favorites and not heavy-favorites 
## median.maxodd <- dt[,max(winodds), by=raceid][,median(V1)]
## ids.small <- dt[,max(winodds), by=raceid][V1 <= median.maxodd, raceid]
## ids.large <- dt[,max(winodds), by=raceid][V1 >= median.maxodd, raceid]
## dt.small <- dt[raceid %in% ids.small]
## dt.large <- dt[raceid %in% ids.large]
## out.small <- mle(minuslogl=function(delta, t) -loglike.sal(delta, t, dt.small),
##             start=list(delta=0.5, t=-0.02) , method="L-BFGS-B",
##             lower=c(0, -2), upper=c(2, 2))
## out.large <- mle(minuslogl=function(delta, t) -loglike.sal(delta, t, dt.large),
##                  start=list(delta=0.5, t=-0.02) , method="L-BFGS-B",
##                  lower=c(0, -2), upper=c(2, 2))

## ##dt <- dt[raceid in 1:1000] 
## ##system.time(
## t0 <- date()
## out <- mle(minuslogl=function(delta, t) -loglike.sal(delta, t, dt),
##             start=list(delta=0.5, t=-0.02) , method="L-BFGS-B",
##             lower=c(0, -2), upper=c(2, 2))
## t1 = date()

## foo <- mcparallel(
##     out <- mle(minuslogl=function(a,t) -loglike.rdu2(a, t, dt),
##                start=list(a=0.5, t=-0.02) , method="L-BFGS-B",
##                lower=c(0.001, -5), upper=c(5, 5))
## )

## bar <- mccollect(foo)[[1]]

## rand <- rbinom(1, 1000, 0.5)
## typ <- "out-rdu-prelec-lbfg"
## ##typ <- "out-rdu-pow-lbfg"
## ##fl <- paste("logs/", typ, rand, ".dat", sep="")
## fl <- paste(typ, rand, ".dat", sep="")
## ##fl.txt <- paste("logs/", typ, rand, ".txt", sep="")
## out = list("t0" = t0, "t1" = t1, "out" = out)
## save(out, file=fl)
## ##dput(summary(out), file=fl.txt)

### vyong test

## no need to run this function again, as dt.vuong already exists in logs,
## it won't run if out.eu, out.r and out.s are not loaded in memory  
save.dt.vuong.first.time <- function() {
    dt.v <- as.data.table(as.data.frame(df.salanie))
    eu.par <- out.eu@coef # t 
    rdu.par <- out.r@coef # a, b, t
    sal.par <- out.s@coef #delta, t
    ## eu
    t <- eu.par[1]
    dt.v[ ,utils:=cara_cpp(winodds+1, t)]
    dt.v[, ps.eu:=solve.ps.eu2(utils), by=raceid]
    dt.v[winner==1, sum(log(ps.eu))] * (-2) ## -2 log L: 141080.8 
    ## rdu 
    a <- rdu.par[1]
    b <- rdu.par[2]
    t <- rdu.par[3]
    dt.v[ ,utils:=cara_cpp(winodds+1, t)]
    dt.v[, ps.rdu:=solve.ps.rdu2(a, b, utils), by=raceid]
    dt.v[winner==1, sum(log(ps.rdu))] * (-2) ## -2 log L: 141073.1 
    ## saliance
    delta <- sal.par[1]
    t <- sal.par[2]
    dt.v[ ,utils:=cara_cpp(winodds+1, t)]
    dt.v[, ps.sal:=solve.ps.sal(delta, utils), by=raceid]
    dt.v[winner==1, sum(log(ps.sal))] * (-2) ## -2 log L: 140269.9 
    dt.vyong <- dt.v
    save(dt.vyong, file="logs/data-table-4-vyongtest.dat")
}

vuong.test <- function() {
    ##dt.vyong already saved
    load("logs/data-table-4-vyongtest.dat")
    dt.v <- dt.vyong
    dt.v[winner==1, sum(log(ps.eu))] * (-2) ## -2 log L: 141080.8 
    dt.v[winner==1, sum(log(ps.rdu))] * (-2) ## -2 log L: 141073.1 
    dt.v[winner==1, sum(log(ps.sal))] * (-2) ## -2 log L: 140269.9 
    ## eu vs rdu
    ## z = LR / (sqrt(n) * omega)
    ## LR = L1 - L2 - (K1 - K2) * log(n)/2
    ## li = log (f1/f2) and omega^2 is var(li^2) 
    ll.eu <- dt.v[winner==1, sum(log(ps.eu))]
    ll.rdu <- dt.v[winner==1, sum(log(ps.rdu))]
    ll.sal <- dt.v[winner==1, sum(log(ps.sal))]
    omega.sal.rdu <- dt.v[winner==1, sqrt(var(log(ps.sal/ps.rdu)))]
    n <- dt.v[,length(unique(raceid))]
    z.sal.rdu <- (ll.sal - ll.rdu) / (sqrt(n) * omega.sal.rdu) ## 15.85127
}


for.crossvalidation <- function() {
    dt <- as.data.table(as.data.frame(df.salanie))
    outs <- list()
    for (yr in dt[,unique(year)]) {
        print(yr)
        dt.sub <- dt[year != yr, ]
        outs[yr] <- mle(minuslogl=function(delta, t) -loglike.sal(delta, t, dt.sub),
                        start=list(delta=0.5, t=-0.02) , method="L-BFGS-B",
                        lower=c(0, -2), upper=c(2, 2))
        print(outs[yr])
    }
    outs
}


## ========================= playground

## t=0.4
## a <- 1.2 
## dt[ ,utils:=cara(winodds+1, t)]
## dt[raceid==999, ps:=solve.ps(a, utils), by=raceid]

## xs <- 1:10

## outs <- list()
## for (yr in dt[,unique(year)]) {
##     print(yr)
##     dt.sub <- dt[year==yr, ]
##     outs[yr] <- mle(minuslogl=function(delta, t) -loglike.sal(delta, t, dt.sub),
##                       start=list(delta=0.5, t=-0.02) , method="L-BFGS-B",
##                       lower=c(0, -2), upper=c(2, 2))
##     print(outs[yr])
## }

## ds <- c()
## ts <- c()
## for (yr in 86:95) {
##     ds <- c(ds, unlist(outs[[yr]]@coef[1]))
##     ts <- c(ts, outs[[yr]]@coef[2])
## }


## load("out-rdu-salanie487.dat")
## out.r <- out$out



