source("util.R")

library(data.table)
library(profvis)

options( stringsAsFactors=F)

dir <-  "~/Documents/data/prediction_data/betfair/horses_2013/" ## office computer
files <- system(paste("ls", dir), intern=T)
          
i <- 16 
tm <- proc.time()

df <- data.table(read.csv(paste(dir, files[i], sep="")))

setnames(df, names(df), sapply(names(df), lowercase))

df[,":="(settled_date = dmy_hms(settled_date),
         scheduled_off = dmy_hms(scheduled_off),
         actual_off = dmy_hms(actual_off),
         latest_taken = dmy_hms(latest_taken),
         first_taken = dmy_hms(first_taken))]

df <- df[in_play != "IP"] # only gambles before the start of the game 

df <- df[odds > 1] # some horsese have odds exactly 1 which does not make sense.  

df[, ":="(p.mean = mean(1/odds),
          p.wght = .SD[, sum(volume_matched/odds)/sum(volume_matched)], 
          p.first = 1/.SD[which.min(first_taken), odds],
          p.last = 1/.SD[which.max(latest_taken), odds]),
          by = c('event_id', 'selection_id')]


## returns a list : [(x1, y1),.., (x_nbin, y_nbin)],
## x_i is the number of horses whose subjective probability of wining is between (i-1)/nbin and i/nbin
## and y_i is number of times horses from this group won.
get.bins <- function(df, nbin=10, prob='p.mean') {
    get.bins.helper(df[, .(p=get(prob)[1], win_flag=win_flag[1]),
                        by=c('event_id','selection_id')],
                        nbin)
}

## df: data.table with 2 columns named p and win_flag
get.bins.helper <- function(df, nbin=10) {
    aux <- function(inds, win_flag) {
        for (ind in inds)
            bins[[ind]] <- bins[[ind]] + c(1, win_flag) 
        assign("bins", bins, envir=parent.frame())
    }

    bins <- list()
    for(i in 1:nbin) {
        bins[[i]] <- c(0, 0)
    }

    aux(sapply(df[win_flag==1, p], which.bin, nbin=nbin),
        win_flag=1)

    aux(sapply(df[win_flag==0, p], which.bin, nbin=nbin),
        win_flag=0)

    return (bins)
}

which.bin <-  function(p, nbin=10) {
   seq(1, nbin)[min(floor(nbin*p) + 1, nbin)] 
}

emperical.probs <- function(df, nbin=10, prob="p.mean") {
    bins <- as.data.table(t(as.data.table(get.bins(df, nbin, prob))))
    bins[,":=" (ind=1:.N,  p.obj=V2/V1)]
    bins[,p.sub := (2*ind-1)/(2*.N)]
    return (bins)
}

plt <- function(df, nbin) {
    par(mfrow=c(2, 2))
    for (p in c("p.mean", "p.wght", "p.first", "p.last")) {
        bins <- emperical.probs(df, nbin, p)
        with(bins, plot(p.obj ~ p.sub, xlim=c(0, 1), ylim=c(0, 1), main=p, axes=F))
        abline(a=0, b=1)
        abline(h=1)
        abline(v=1)
        abline(h=0)
        abline(v=0)
    }
    par(mfrow=c(1, 1))
    return (NULL)
}

plt(df, 20)
print(length(df[, unique(event_id)]))
time.taken <- proc.time() - tm
print(time.taken)

## playground


## p <- seq(0,1, 0.01)
## q <- sapply(p, fun, sig=0.01)
## plot(c(0.5), xlim = c(0, 1), ylim=c(0, 1))
## lines(q~p, type="l")
## lines(p^1~p, col="red")
## abline(h=0.5)
## abline(v=0.5)

pdf("~/work/misc/conferences-workshops-talks/longshotbias/resources/betfair.pdf")
par(mfrow=c(1, 1))
p <- "p.last"
nbin <- 25 
bins <- emperical.probs(df, nbin, p)
u <- 0.5 
with(bins, plot(p.sub ~ p.obj, xlim=c(0, u), ylim=c(0, u), axes=F, xlab="p", ylab=expression(hat(p))) )
axis(1, at=seq(0, 0.5, 0.1))
axis(2, at=seq(0, 0.5, 0.1))
abline(a=0, b=1)
dev.off()
