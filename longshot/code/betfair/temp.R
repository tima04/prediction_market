library(data.table)
library(lubridate)

options( stringsAsFactors=F)
df <- data.table(read.csv("~/Documents/prediction_data/betfair/horses_2014/140106.csv"))

nms <- sapply(names(df), lowercase)
names(nms) <- NULL
names(df) <- nms

## dt <- df[,11:18,with=F][,-c(5, 6), with=F]
## dt <- dt[1:10000]


## foo <- df[isthere("wadswick", selection)] 

## utility functions

## isthere("ab", c("Aab", "daf")) == c(TRUE, FALSE)
isthere <- function(pat, xs) {
    match <- function(x)
        length(grep(pat, x, ignore.case = T, value = T)) > 0
    rslt <- unlist(Map(match, xs))
    names(rslt) <- NULL
    return (rslt)
}

## lowercase("Ab_CdE") ==  "ab_cde"
lowercase <- function(chars) {
    big <- "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    small <- "abcdefghijklmnopqrstuvwxyz"
    rslt <- ""
    for (char in strsplit(chars, "")[[1]]) {
            pos <- regexpr(char, big)
        if (pos %in% seq(1,26)) 
            rslt <- paste(rslt, (strsplit(small, "")[[1]][pos]), sep="")
        else 
            rslt <- paste(rslt, char, sep="")
    }
    rslt
}

to.seconds <- function(times) {
   as.numeric(dmy_hms(times)) 
}


## playground
ids <- df[1:100, unique(event_id)]

cat("\n\n\n")
n <- 25
foo <- df[event_id==ids[n]] 
setkey(foo, first_taken)
##foo <- df[event_id==ids[n] & in_play=="PE"]
##bar <- foo[, c(10, 11, 12, 13, 14, 17, 18), with=F]
foo[, list(mean=mean(odds), median=median(odds), sd=sd(odds)), by=selection]
print(paste("winner is:",  bar[win_flag==1, unique(selection)],
            "prob:", sum(bar[, 1/mean(odds), by=selection][,V1])))

horse <- "amplefo"
baz <- foo[isthere(horse, selection)]
ind <- 1 + sum(baz[,in_play] == "PE")
plot(100/baz[,odds], type="l", log='y')
abline(v=ind)


#get all races of the horse wolfs Girl
horse = '"wolfs girl"'
dir <- "~/Documents/prediction_data/betfair/horses_2014/*csv"
cmd <- paste("grep -ih", horse, dir, "> ", "horse.csv")
system(cmd)
horse <- data.table(read.csv("horse.csv", header = F)) 
names(horse) <- names(df)

## binning 
bar <- bar[1:50]
nbin <- 10

bin.count <- function(df, nbin=10) {
    df[, ":="(p = mean(1/odds)), by=selection_id]
    bin.count.helper(df[, list(p=p[1], win=win_flag[1]), by=selection_id][, 2:3, with=F],
              nbin)
}

## df: data.table with 2 columns named p and win
bin.count.helper <- function(df, nbin=10) {
    aux <- function(inds, is.winner) {
        for (ind in inds)
            bins[[ind]] <- bins[[ind]] + c(1, is.winner) 
        assign("bins", bins, envir=parent.frame())
    }

    bins <- list()
    for(i in 1:nbin) {
        bins[[i]] <- c(0, 0)
    }

    aux(sapply(df[win==1, p], which.bin, nbin=nbin),
        is.winner=1)

    aux(sapply(df[win==0, p], which.bin, nbin=nbin),
        is.winner=0)

    return (bins)
}

which.bin <-  function(p, nbin=10) {
   seq(1, nbin)[min(floor(nbin*p) + 1, nbin)] 
}

emperical.probs <- function(df, nbin=10) {
    bins <- as.data.table(t(as.data.table(bin.count(df, nbin))))
    bins[,":=" (ind=1:.N,  p.obj=V2/V1)]
    bins[,p.sub := (2*ind-1)/(2*.N)]
    return (bins)
}

plt <- function(df, nbin) {
    bins <- emperical.probs(df, 10)
    with(bins, plot(p.sub ~ p.obj, xlim=c(0, 1), ylim=c(0, 1)))
    abline(a=0, b=1)
}
plt(df, 20)

#############################

## unbiased but random beliefs 2 horses
fun <- function(p = 0.4, n=1000, sig=0.1) {
    n1 <- c() 
    n2 <- c() 
    for (i in 1:n) {
       q <- rnorm(1, p, sig)
       if (ev(q, sum(n1), sum(n2)) == 1)
           n1 <- c(n1, 1)
       else ## d == 2 
           n2 <-  c(n2, 1)
    }
    ##print(c(sum(n1), sum(n2)))
    p.implied <- sum(n1) / (sum(n1) + sum(n2))
    p.implied
}

ev <- function(p, n1, n2) {
    if ((n1+n2) == 0)
        return (rbinom(1, 1, 0.5))  
    if (p > n1/(n1+n2))
        return (1)
    return (2)
}

plt <- function(sig=0.1, n=10000, length=100) {
    ps <- seq(0.1, 1, length=length) 
    qs <- sapply(ps, fun, sig=sig) 
    plot(qs~ps, xlim=c(0, 1), ylim=c(0, 1), type="p", col='red')
    ## points(ps, ps, col="red")
    abline(b=1)
    abline(a=0, b=1)
    abline(h=0.5)
    abline(v=0.5)
}
plt(sig=100, n=10000, length=20)

#end ###########################

#############################

## unbiased but random beliefs nh horses
fun <- function(nh=4, n=1000, sig=0.1) {
    ps <- c(0.1, 0.2, 0.3, 0.4)
    hs <- rep(0, nh)  
    for (i in 1:n) {
       q <- rnorm(1, p, sig)
       if (ev(q, sum(n1), sum(n2)) == 1)
           n1 <- c(n1, 1)
       else ## d == 2 
           n2 <-  c(n2, 1)
    }
    ##print(c(sum(n1), sum(n2)))
    p.implied <- sum(n1) / (sum(n1) + sum(n2))
    p.implied
}

ev <- function(p, n1, n2) {
    if ((n1+n2) == 0)
        return (rbinom(1, 1, 0.5))  
    if (p > n1/(n1+n2))
        return (1)
    return (2)
}

plt <- function(sig=0.1, n=10000, length=100) {
    ps <- seq(0.1, 1, length=length) 
    qs <- sapply(ps, fun, sig=sig) 
    plot(qs~ps, xlim=c(0, 1), ylim=c(0, 1), type="p", col='red')
    ## points(ps, ps, col="red")
    abline(b=1)
    abline(a=0, b=1)
    abline(h=0.5)
    abline(v=0.5)
}
plt(sig=0.2, n=10000, length=20)

#end ###########################

w <- function(p, a=0.6) {
    exp(-(-log(p))^a)
}

fun <- function(p=0.4, n1=4, n2=6, a=1, n=1000, pt=F) {
    pis <- c()
    for (i in 1:n) {
        pi <- n1/(n1+n2)
        pis <- c(pis, pi)
        choice = ifelse(pt, pt(pi, p, a), eu(pi, p))
        if (choice == 1)
            n1 <- n1+1
        else
            n2 <- n2+1
    }
    pis[n]
}

pt <- function(pi, p, a) {
    ifelse(w(p, a)/pi > w(1-p, a)/(1-pi), 1, 2)
}

eu <- function(pi, p, a=1, b=0, c=0) {
    ifelse(p*u(1/pi, p, a, b, c) > p*u(1/(1-pi), 1-p, a, b, c), 1, 2)
}

u <- function(x, p, a=1, b=0, c=0) {
    p*(a*x + b*x^2 + c*x^3)
}


