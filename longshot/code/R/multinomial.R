library(data.table)
library(nleqslv)
data.dir <- "../data/gandhi_restd/" 

df.clm <- data.table(read.csv(paste(data.dir, "dataset_clm.csv", sep=""))) ## maiden horses
df.mcl <- data.table(read.csv(paste(data.dir, "dataset_mcl.csv", sep=""))) ## non-maiden(experienced) 

## combining 
df.clm[, maiden:=1]
df.mcl[, maiden:=0]
df <- as.data.table(rbind(df.clm, df.mcl))
setkey(df, raceid, winodds)

## EU with Cara
## let phi = theta * a and gi = 1/(exp(phi) - exp(-phi * Ri))
## then pi = gi/(sum(gj)) 
## loglike = sum(log(pk)) where pk is the sucesssfull horse in i'th group.

loglike <- function(phi, dt) {
   dt[, gi:= 1/(exp(phi) - exp(-phi*winodds))] 
   dt[, hi := gi/sum(gi), by=raceid]
   return (dt[winner==1, sum(log(hi))])
}

dt <- as.data.table(as.data.frame(df[maiden==0]))
#loglike(1, dt)
baz <- optimise(f=function(x) -loglike(x, dt), interval = c(-1, 1))

## pt with cara and power weighting parameters.

## let in a race p_i: the probability of i'th horse wining & n: num of horses
## w = p_i^a (1 - exp(-t*R_i)) + (1 - p_i)^b (1 - exp(-t)) (pars: a, b, t) gives n equation
## and sum(ps) = 1, gives 1 equation, these can be used to solve for ps.
solve.ps <- function(a, b, t, Rs, epsilon=10e-9) {
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

## phi: (alpha, beta, theta)
loglike.pt <- function(phi, dt) { 
    a <- phi[1]
    b <- phi[2]
    t <- phi[3]
    dt[, ps:=solve.ps(a, b, t, winodds), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}


## method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
numrace <- 500 
ids <- df[,unique(raceid)]
dt <- as.data.table(as.data.frame(df[maiden==0]))
#dt <- dt[raceid %in% ids[1:numrace]]
t0 <- Sys.time()
par.jul <- c(-0.72, 1.162, 0.318)
out <- optim(par=c(0.5, 0.5, 0.5), fn=function(x) -loglike.pt(x, dt), gr=NULL, method="L-BFGS-B")
#out <- optim(par=c(0.5, 0.5, 0.5), fn=function(x) -loglike.pt(x, dt), gr=NULL, method="Nelder-Mead")
t1 <- Sys.time()
#out
tm <- t1 - t0
fl <- paste("out", rbinom(1, 100, 0.5), ".txt", sep="")
#save(tm, out, file=fl)
writeLines(paste(out, paste("time: ", tm), sep="\n"), fl)
#writeLines(paste("time: ", tm, sep=""), fl)

xs <- seq(-0.05, 0.05, .001)
ys <- sapply(xs, loglike, dt=dt)
##plot(ys~xs, type="l")

## playground
foo <- function(x, theta) {
    (0.5 * (x + 100)^theta + 0.5 * x^theta)^(1/theta) - x
} 

bar <- function(xs) {
    return (c(sum(xs^2) - 1, xs[1] -0.5, xs[2] + xs[3]))
}
## rslt <- nleqslv(c(10, 100, 10), bar)$x
## rslt
## sum(rslt^2)


## {% Formulate a variation of EU where both regret and disappointment are incorporated, and show how particular assumptions on the form of utility lead to empirical predictions such as the Allais paradox. %}
## Laciana, Carlos E. & Elke U. Weber (2008) “Correcting Expected Utility for Comparisons between Alternative Outcomes: A Unified Parameterization of Regret and Disappointment,” Journal of Risk and Uncertainty 36, 1–17
