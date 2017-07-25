library(data.table)
library(stats4)
library(nleqslv)
data.dir <- "../data/gandhi_restd/" 

df.clm <- data.table(read.csv(paste(data.dir, "dataset_clm.csv", sep=""))) ## maiden horses
df.mcl <- data.table(read.csv(paste(data.dir, "dataset_mcl.csv", sep=""))) ## non-maiden(experienced) 

## combining 
df.clm[, maiden:=1]
df.mcl[, maiden:=0]
df <- as.data.table(rbind(df.clm, df.mcl))
setkey(df, raceid, winodds)

utility <- function(xs, theta) {
    (1 - exp(-theta*xs))/theta
}

regret <- function(xs, reg.par) {
   1 - reg.par^xs 
}

## sum.other(c(x1, x2, x3)) = c(x2+x3, x1+x3, x2+x3)
sum.other <- function(xs) {
    sum(xs) - xs
}

solve.ps.reg <- function(reg.par, theta, Rs, epsilon=10e-9) {
        n <- length(Rs)

        fun <- function(ps.plus.w) { # ps.plus.w = c(ps, w)
            ps <- ps.plus.w[1:n]
            w = ps.plus.w[n+1]
            u <- utility
            ## ws is (n-1)[u(-1)(1-p_i) + u(R_i) + regret(u(Ri) - u(-1))] +
            ## sum(pj * regret(u(-1) - u(R_j))) where j in {1,..,n} - {i}
            u.rs <- utility(Rs, theta)
            u.min.1 <- utility(-1, theta)
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

loglike.reg <- function(phi, dt) {
    reg.par <- phi[1]
    theta <- phi[2]
    dt[, ps:=solve.ps.reg(reg.par, theta, winodds), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

loglike.reg2 <- function(reg.par, theta, data=dt) {
    data[, ps:=solve.ps.reg(reg.par, theta, winodds), by=raceid]
    rslt <- dt[winner==1, sum(log(ps))]
    rslt
}

## method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
dt <- as.data.table(as.data.frame(df))
## numrace <- 500 
## ids <- df[,unique(raceid)]
## n <- as.integer(0.01 * length(ids))
## sids <- sample(ids, n)
## dt <- dt[raceid %in% sids]
out <- mle(minuslogl=function(r, t) -loglike.reg2(r, t), start=list(r=0.5, t=0.5),
           method="Nelder-Mead")
rand <- rbinom(1, 1000, 0.5)
typ <- "out-full-like-reg-"
fl <- paste("logs/", typ, rand, ".dat", sep="")
fl.txt <- paste("logs/", typ, rand, ".txt", sep="")
save(out, file=fl)
dput(summary(out), file=fl.txt)
