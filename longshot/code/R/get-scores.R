library(data.table)
library(stats4)
library(nleqslv)

source("util.R")
source("pt.R")
source("regret2.R")
source("salience2.R")

data.dir <- "../data/gandhi_restd/" 
df <- fread(paste(data.dir, "dataset_full.csv", sep=""))
setkey(df, raceid, winodds)

dt <- as.data.table(as.data.frame(df))

get.pt.score <- function(pars, winodds, winner, prelec=FALSE) {
    a <- pars[1]
    b <- pars[2]
    t <- pars[3]
    #browser()
    if (prelec) 
        solve = solve.ps2
    else
        solve = solve.ps
    #solve <- ifelse(prelec, solve.ps2, solve.ps)
    ps <- solve(a, b, t, winodds)
    return (log(sum(ps*winner)))
}

get.reg.score <- function(pars, winodds, winner, prelec=FALSE) {
    reg.par <- pars[1]
    t <- pars[2]
    ps <- solve.ps.reg(reg.par, t, winodds)
    return (log(sum(ps*winner)))
}

get.sal.score <- function(pars, winodds, winner, prelec=FALSE) {
    delta <- pars[1]
    t <- pars[2]
    ps <- solve.ps.st(delta, t, winodds)
    return (log(sum(ps*winner)))
}

load("logs/out-full-like-reg-505.dat")
out.reg <- out
load("logs/out-full-like-sal-506.dat")
out.sal<- out
load("out-full2-like-524.dat")
out.pt<- out
load("logs/out-pt2-power-480.dat")
out.pt.pow <- out
load("logs/out-pt2-prelec-531.dat")
out.pt.prl <- out
load("logs/out-pt2-prelec-lbfg-502.dat")
out.pt.prl.lbfg <- out

dt[, ":=" (pt.pw.scr = get.pt.score(out.pt.pow@coef, winodds, winner), 
           pt.prl.scr = get.pt.score(out.pt.prl@coef, winodds, winner),
           reg.scr = get.reg.score(out.reg@coef, winodds, winner),
           sal.scr = get.sal.score(out.sal$par, winodds, winner)), 
           by=raceid]


######### playground
dt2 <- dt[pt.pw.scr > -Inf] 
fun <- median
colSums(dt2[, .(pt=pt.prl.scr[1], reg=reg.scr[1], sal=sal.scr[1]), by=raceid])

foo <- dt2[pt.pw.scr < reg.scr]
ids <- foo[,unique(raceid)]

bar <- dt2[pt.pw.scr > reg.scr]
ids2 <- bar[,unique(raceid)]

n <- 7 
foo[raceid==ids[n]]
bar[raceid==ids2[n]]

dt3 <- dt[pt.pw.scr == -Inf]


ids <- dt3[,unique(raceid)]

n <-2  
dt3[raceid == ids[n]]

