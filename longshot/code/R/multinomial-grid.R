library(data.table)
library(nleqslv)
library(stats4)

data.dir <- "../data/gandhi_restd/" 

df.clm <- data.table(read.csv(paste(data.dir, "dataset_clm.csv", sep=""))) ## maiden horses
df.mcl <- data.table(read.csv(paste(data.dir, "dataset_mcl.csv", sep=""))) ## non-maiden(experienced) 
df.full <- fread("../data/gandhi_restd/dataset_full.csv") ## somehow contains more then clm + mcl

## combining 
df.clm[, maiden:=1]
df.mcl[, maiden:=0]
df <- as.data.table(rbind(df.clm, df.mcl))
setkey(df, raceid, winodds)


utility <- function(xs, theta) {
    (1 - exp(-theta*xs))/theta
}

w.prelec <- function(ps, alpha) {
    exp(-(-log(ps))^alpha)
}

## pt with cara and power weighting parameters.

## let in a race p_i: the probability of i'th horse wining & n: num of horses
## w = p_i^a (1 - exp(-t*R_i)) + (1 - p_i)^b (1 - exp(-t)) (pars: a, b, t) gives n equation
## and sum(ps) = 1, gives 1 equation, these can be used to solve for ps.
solve.ps <- function(a, b, t, Rs, epsilon=10e-100) {
    n <- length(Rs)
    fun <- function(xs) { # xs = c(ps, w)
        ps <- xs[1:n]
        w = xs[n+1]
        ## ws is ps^a * (1 - exp(-t*Rs)) + (1 - ps)^b * (1 - exp(-t))
        rslt <- ps^a * (1 - exp(-t*Rs))/t + (1 - ps)^b * (1 - exp(t))/t - w
        return (c(rslt, sum(ps)-1))
    }
    ps <- rep(1/n, n) # initial value
    w <- 0.5
    rslt <- nleqslv(c(ps, w), fun)$x[1:n]
    rslt <- ifelse(rslt <= 0, epsilon, rslt)
    rslt
}

## prelec weighting function
solve.ps2 <- function(a, b, t, Rs, epsilon=10e-10) {
    n <- length(Rs)
    fun <- function(xs) { # xs = c(ps, w)
        ps <- xs[1:n]
        w = xs[n+1]
        ## ws is w.prelec(ps, a) * utility(Rs, t) + w.prelec((1 - ps),b) * utility(-1, t) 
        rslt <-  w.prelec(ps, a) * utility(Rs, t) + w.prelec((1 - ps),b) * utility(-1, t) - w
        return (c(rslt, sum(ps)-1))
    }
    ps <- rep(1/n, n) # initial value
    w <- 0.5
    rslt <- nleqslv(c(ps, w), fun)$x[1:n]
    rslt <- ifelse(rslt <= 0, epsilon, rslt)
    rslt
}

## phi: (alpha, beta, theta)
loglike.pt <- function(phi, dt) { 
    a <- phi[1]
    b <- phi[2]
    t <- phi[3]
    dt[, ps:=solve.ps(a, b, t, winodds), by=raceid]
    ##dt[, ps:=solve.ps2(a, b, t, winodds), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

loglike.pt2 <- function(a, b, t, dt) { 
    ##dt[, ps:=solve.ps(a, b, t, winodds), by=raceid]
    dt[, ps:=solve.ps2(a, b, t, winodds), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

#### run pt on high and low longshot (cutoff is median) 
## dt <- as.data.table(as.data.frame(df[maiden==1]))
##dt <- as.data.table(as.data.frame(df))
dt <- as.data.table(as.data.frame(df.full))
## ids <- dt[,unique(raceid)]
## n <- as.integer(0.01 * length(ids))
## sids <- sample(ids, n)
## dt <- dt[raceid %in% sids]
## median.max.odd <- dt[, .(odd = max(winodds)), by=raceid][,median(odd)]
## id.h <- dt[, .(odd = max(winodds)), by=raceid][odd > median.max.odd, raceid]
## id.l <- dt[, .(odd = max(winodds)), by=raceid][odd < median.max.odd, raceid]
## dt.h <- dt[raceid %in% id.h]
## dt.l <- dt[raceid %in% id.l]
## t0 <- Sys.time()
## out.h <- optim(par=c(0.5, 0.5, 0.5), fn=function(x) -loglike.pt(x, dt.h), gr=NULL, method="Nelder-Mead")
## out.l <- optim(par=c(0.5, 0.5, 0.5), fn=function(x) -loglike.pt(x, dt.l), gr=NULL, method="Nelder-Mead")
out <- mle(minuslogl=function(a,b,t) -loglike.pt2(a, b, t, dt), start=list(a=0.5, b=0.5, t=0.5),
           ,method="Nelder-Mead")
## out.h <- mle(minuslogl=function(a,b,t) -loglike.pt2(a, b, t, dt.h), start=list(a=0.5, b=0.5, t=0.5),
##            ,method="Nelder-Mead")
## out.l <- mle(minuslogl=function(a,b,t) -loglike.pt2(a, b, t, dt.l), start=list(a=0.5, b=0.5, t=0.5),
##            ,method="Nelder-Mead")


rand <- rbinom(1, 1000, 0.5)
typ <- "out-full-like-prelec-"
fl <- paste("logs/", typ, rand, ".dat", sep="")
fl.txt <- paste("logs/", typ, rand, ".txt", sep="")
save(out, file=fl)
dput(summary(out), file=fl.txt)
#save(out, out.h, out.l, file=fl)
## fl.l <- paste(typ, "-l-", rand, ".txt", sep="")
## fl.h <- paste(typ, "-h-", rand, ".txt", sep="")

## writeLines(paste(out.l), fl.l)
## writeLines(paste(out.h), fl.h)
## writeLines(paste(out.l), fl.l)


#####                            playground
##=======================================================================

## util <- function(x, t) {
##     (1 - exp(-t*x))/t
## }

## w.p <- function(p, a) {
##     p^a
## }

## w.n <- function(p, b) {
##     p^b
## }

## pt.val <- function(xs, ps, t, a, b) {
##     us <- util(xs, t)
##     w.p(ps[1], a) * us[1] + w.n(ps[2], b) * us[2]
## }

## xs <- c(10^3, -1)
## p <- 10^(-3) 
## ps <- c(p, 1-p)
## pt.val(xs, ps, -0.033, 1.26, .137)
## pt.val(xs, ps, -0.047, 1.25, .27)

## x1 <- 1.3 
## x2 <- 5.3 
## util(x1, -0.033)
## util(x2, -0.047)

## t <- -0.072
## a <- 1.162
## b <- 0.318
## n <- 3 
## xs <- c(10^n, -1)
## p <- 10^(-n)
## ps <- c(p, 1-p)
## pt.val(xs, ps, t, a, b)
## util(3, t)

## Estimate  Std. Error
## a  1.26305545 0.052771350
## b  0.13747371 0.073848764
## t -0.03266529 0.004557875

## -2 log L: 51731.43 
## > summary(out.l)
## Maximum likelihood estimation

## Call:
## mle(minuslogl = function(a, b, t) -loglike.pt2(a, b, t, dt.l), 
##     start = list(a = 0.5, b = 0.5, t = 0.5), method = "Nelder-Mead")

## Coefficients:
##      Estimate Std. Error
## a  1.24607074 0.08661108
## b  0.26999406 0.13655491
## t -0.04756173 0.01077767

