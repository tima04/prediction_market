library(data.table)
library(stats4)
library(nleqslv)
library(Rcpp)
library(parallel)

source("util.R")
sourceCpp("utils.cpp") 

regret <- function(xs, reg.par) {
   1 - reg.par^xs 
}

## sum.other(c(x1, x2, x3)) = c(x2+x3, x1+x3, x2+x3)
sum.other <- function(xs) {
    sum(xs) - xs
}

solve.ps <- function(reg.par, utils, epsilon=10e-9) {
        n <- length(utils)
        fun <- function(xs) { # xs = c(ps, w) 
            ps <- xs[1:n]
            w <- xs[n+1]
            ## ws is (n-1)[u(-1)(1-p_i) + u(R_i) + regret(u(Ri+1) - u(-1))] +
            ## sum(pj * regret(u(-1) - u(R_j+1))) where j in {1,..,n} - {i}
            ## since u(-1) = 0, in new parameteraziation this is equivalent to 
            ## ws[i] is (n-1)[ pi * util[i] + regret(util[i]) +
            ## sum(pj * regret(-util[j])) where j in {1,..,n} - {i}
            ## rslt <- (n-1) * (utils + regret(utils, reg.par)) +
            ##     sum.other(ps * regret(-utils, reg.par)) -
            ##     w
            rslt <- (ps * utils + regret(utils, reg.par)) +
                sum.other(ps * regret(-utils, reg.par)) -
                w
            return (c(rslt, sum(ps)-1))
        }
        ps <- rep(1/n, n) # initial value
        w <- 0.5
        rslt <- nleqslv(c(ps, w), fun)$x[1:n]
        #rslt <- dfsane(c(ps, w), fun, control=list(trace=FALSE))$par[1:n]
        rslt
}

solve.ps.old <- function(reg.par, theta, Rs, epsilon=10e-9) {
        n <- length(Rs)
        fun <- function(ps.plus.w) { # ps.plus.w = c(ps, w)
            ps <- ps.plus.w[1:n]
            w = ps.plus.w[n+1]
            u <- cara 
            ## ws is (n-1)[u(-1)(1-p_i) + u(R_i) + regret(u(Ri+1) - u(-1))] +
            ## sum(pj * regret(u(-1) - u(R_j+1))) where j in {1,..,n} - {i}
            u.rs <- cara(Rs+1, theta)
            u.min.1 <- cara(-1, theta)
            regs1 <- regret(u.rs - u.min.1, reg.par)
            regs2 <- regret(u.min.1 - u.rs, reg.par)
            rslt <- (n-1) * (u.min.1 * (1 - ps) + u.rs + regs1) +
                sum.other(ps * regs2) - w
            return (c(rslt, sum(ps)-1))
        }
        ps <- rep(1/n, n) # initial value
        w <- 0.5
        rslt <- nleqslv(c(ps, w), fun)$x[1:n]
        rslt <- ifelse(rslt <= 0, epsilon, rslt)
        rslt
}

loglike.reg <- function(reg.par, theta, dt) {
    dt[ ,utils:=cara_cpp(winodds+1, t)]
    dt[, ps:=solve.ps(reg.par, utils), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

dt <- as.data.table(as.data.frame(df))
out <- mle(minuslogl=function(r, t) -loglike.reg(r, t), start=list(r=0.5, t=0.5),
        method="Nelder-Mead")
rand <- rbinom(1, 1000, 0.5)
typ <- "out-reg2-fltered-wincorrected-"
fl <- paste("logs/", typ, rand, ".dat", sep="")
fl.txt <- paste("logs/", typ, rand, ".txt", sep="")
save(out, file=fl)
dput(summary(out), file=fl.txt)

#### playground
require("BB")
froth <- function(p){
    r <- rep(NA, length(p))
    r[1] <- -13 + p[1] + (p[2] * (5 - p[2]) - 2)  * p[2]
    r[2] <- -29 + p[1] + (p[2] * (1 + p[2]) - 14) * p[2]
    r
}

fun <- function(xs) {
    r <- rep(NA, length(xs))
    r[1] <- xs[1]^2 - xs[2]^2 
    r[2] <- xs[2] - 3
    r
}

p0 <- rep(0, 2)
dfsane(par = p0, fn = fun)


par[1] -4.990589 -1.448398$residual[1] 10.42073$fn.reduction[1] 17.04337$feval[1] 131$iter[1] 106$convergence

dfsane(par = p0, fn = froth, control = list(trace = FALSE))$par[1] -4.990589 -1.448398$residual[1] 10.42073$fn.reduction[1] 17.04337$feval[1] 131$iter[1] 106$convergence
