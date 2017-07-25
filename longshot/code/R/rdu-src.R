library(data.table)
library(stats4)
library(nleqslv)
source("util.R")


w.prelec <- function(ps, alpha) {
    exp(-(-log(ps))^alpha)
}

## solve.ps <- function(t, Rs, epsilon=10e-10) {
##     n <- length(Rs)
##     us <- cara(Rs + 1, t)
##     u.min.1 <- cara(-1, t)
##     w <- 1/sum(1/(us - u.min.1)) + u.min.1 
##     ps <- (w - u.min.1)/(us - u.min.1)
##     return(ps)
## }

## prelec weighting function
solve.ps <- function(a, t, gamma, Rs, epsilon=10e-10) {
    n <- length(Rs)
    fun <- function(xs) { # xs = c(ps, w)
        ps <- xs[1:n]
        w = xs[n+1]
        as <- a + gamma * seq(0, 1, length=n)
        ## ws is w.prelec(ps, a) * cara(Rs, t) + w.prelec((1 - ps),b) * cara(-1, t) 
        wghts <- w.prelec(ps, as)
        rslt <-  wghts * cara(Rs, t) + (1 - wghts) * cara(-1, t) - w
        return (c(rslt, sum(ps)-1))
    }
    ps <- rep(1/n, n) # initial value
    w <- 0.5
    rslt <- nleqslv(c(ps, w), fun, control = list(allowSingular=TRUE))$x[1:n]
    rslt <- ifelse(rslt <= 0, epsilon, rslt)
    rslt
}

nproblems <- 0 ## number of races not solved
problem.Rs <- list()
solve.ps.catch.error <- function(a, t, gamma, Rs) {
    tryCatch({
        solve.ps(a, t, gamma, Rs)
    }, error = function(e) {
        print(paste(a, t, gamma, Rs, sep=""))
        problem.Rs[nproblems] <- Rs
        nproblems <<- nproblems + 1
        return (NA)
    }
   ) 
}

loglike.rdu <- function(a, t, gamma, dt) { 
dt[, ps:=solve.ps.catch.error(a, t, gamma, winodds), by=raceid]
rslt <- dt[winner==1, sum(log(ps), na.rm=T)]
rslt
}

##dt <- as.data.table(as.data.frame(df))
#out <- mle(minuslogl=function(a,t) -loglike.rdu(a, t, dt), start=list(a=0.5, t=0.5) ,method="L-BFGS-B")
## rand <- rbinom(1, 1000, 0.5)
## typ <- "out-rdu-prelec-lbfg-filter-"
## fl <- paste("logs/", typ, rand, ".dat", sep="")
## fl.txt <- paste("logs/", typ, rand, ".txt", sep="")
## save(out, file=fl)
## dput(summary(out), file=fl.txt)

ids <- df[,unique(raceid)] 
frac <- 0.2
set.seed(10)
sids <- sample(ids, as.integer(frac*length(ids)))
dt <- as.data.table(as.data.frame(df[raceid %in% sids]))
system.time(out <- mle(minuslogl=function(a,t,gamma) -loglike.rdu(a, t, gamma,dt),
                       start=list(a=0.5, t=0.5, gamma=0) ,method="L-BFGS-B"))

## ==========playground
xs <- 1:10
sample(xs, 3)
set.seed(10)
sample(xs, 3)

foo <- function(x) {
    if (x > 0)
        1/x
    else
        stop("error")
}

baz <- 0
bar <- function(x) {
    tryCatch({
        foo(x)
    }, error = function(e) {
       print(e) 
       baz <<- baz + 1
       NA
    }
   ) 
}


result = tryCatch({
    expr
}, warning = function(w) {
    warning-handler-code
}, error = function(e) {
    error-handler-code
}, finally = {
    cleanup-code
}

when rs+1 was used
> summary(out.rdu.src)
Maximum likelihood estimation

Call:
mle(minuslogl = function(a, t, gamma) -loglike.rdu(a, t, gamma, 
    dt), start = list(a = 0.5, t = 0.5, gamma = 0), method = "L-BFGS-B")

Coefficients:
         Estimate  Std. Error
a      0.77123953 0.014681512
t     -0.00147516 0.001273526
gamma  0.13557324 0.018059799

-2 log L: 61313.55 

a <- 0.953
ag <- a + 0.0638
ps <- seq(0, 1, length=100)
plot(w.prelec(ps, a), type="l")
lines(w.prelec(ps, ag), type="l", col="red")

