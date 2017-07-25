library(data.table)
library(stats4)
library(nleqslv)
library(Rcpp)
library(parallel)
library(ggplot2)

source("util.R")
sourceCpp("utils.cpp") 
source("rdu.R")
source("eu.R")
source("salience.R")

df <- as.data.table(as.data.frame(df.salanie))

### extract rdu coeffs
load("outs4cv-rdu.dat")
outs.rdu <- out$out
rdu.coef <- list()
for (yr in 86:95) 
    rdu.coef[yr] <- list(c(outs.rdu[[yr]]@coef))

## sal coeffs
load("outs4cv-salience.dat")
outs.sal <- out.sal4cv
## eu
load("outs4cv-eu.dat")
outs.eu <- out.eu4cv


dt <- as.data.table(as.data.frame(df))
get.cv.rdu <- function() {
    coefs <- rdu.coef
    rslt <- c()
    for (yr in 86:95) {
        rslt <- c(rslt, loglike.rdu(
                            coefs[[yr]]['a'], coefs[[yr]]['beta'], coefs[[yr]]['t'],
                            dt[year==yr]))
    }
    rslt
}

get.cv.eu <- function() {
    rslt <- c()
    for (yr in 86:95) {
        rslt <- c(rslt, loglike.eu(
                            outs.eu[[yr]]@coef, dt[year==yr]))
    }
    rslt
}

get.cv.sal <- function() {
    coefs <- outs.sal 
    rslt <- c()
    for (yr in 86:95) {
        rslt <- c(rslt, loglike.sal(
                            coefs[[yr]]@coef['delta'], coefs[[yr]]@coef['t'], 
                            dt[year==yr]))
    }
    rslt
}

system.time(cv.rdu <- get.cv.rdu())
system.time(cv.eu <- get.cv.eu())
system.time(cv.sal <- get.cv.sal())


## ============================= t test between salience and rdu for each fold
## t test of difference between loglikehood of sal and rdu at each data point 
t.tst.sal.rdu <- function(pval=FALSE) {
    coefs.rdu <- rdu.coef
    coefs.sal <- outs.sal 
    rslt <- c()
    for (yr in 86:95) {
        ## extract rdu coef
        a <- coefs.rdu[[yr]]['a']
        b <- coefs.rdu[[yr]]['beta']
        t.rdu <- coefs.rdu[[yr]]['t']
        ## extract sal coef
        delta <- coefs.sal[[yr]]@coef['delta']
        t.s <- coefs.sal[[yr]]@coef['t']
        rslt <- c(rslt,
                  t.tst.sal.rdu.helper(delta, t.s, a, b, t.rdu, dt[year==yr],
                                       pval=pval))
    }
    rslt
}

t.tst.sal.rdu.helper <- function(delta, t.s, a, b, t.rdu, dt, pval=TRUE) {
    dt[ , utils.rdu:=cara_cpp(winodds+1, t.rdu)]
    dt[ , ps.rdu:=solve.ps.rdu2(a, utils.rdu, beta=b), by=raceid]
    dt[ , utils.sal:=cara_cpp(winodds+1, t.s)]
    dt[ , ps.sal:=solve.ps.sal(delta, utils.sal), by=raceid]
    t.stat <- dt[winner==1, log(ps.sal) - log(ps.rdu)]
    if (pval)
         return (t.test(t.stat, alternative = "greater")$p.val)
    return (sd(t.stat) * sqrt(dt[,.N]))
}

sal.rdu.sd <- t.tst.sal.rdu()

frm <- data.frame(year=86:95, ll=cv.sal-cv.rdu, sd=sal.rdu.sd) 
alpha.sig <- 1 
p <- ggplot(data=frm, aes(x=year, y=ll)) +
    geom_line() +
    geom_point() +
    ##geom_ribbon(aes(ymin=ll-alpha.sig*sd, ymax=ll+alpha.sig*sd), alpha=0.2) +
    theme_minimal() +
    ylim(0, 70) +
    xlab("Year") +
    ylab("LL(salient-theory) - LL(prospect-theory)") +
    ##scale_x_date()
    scale_x_continuous(breaks=seq(86, 95, 1),
                       label=paste("19", seq(86, 95, 1), sep=""))
pdf("../../presentation/images/cv.pdf")
print(p)
dev.off()

## =========================== Playground ==================================

for (y in 86:95)
    print(outs.sal[[y]]@coef)
yr=92
coefs <- outs.sal 
loglike.sal(coefs[[yr]]@coef['delta'], coefs[[yr]]@coef['t'], dt[year==yr])
