library(data.table)

data.dir <- "../../data/" 
df <- fread(paste(data.dir, "dataset_full.csv", sep=""))
setkey(df, raceid, winodds)

## remove races for which no winner in the data eg: raceid 19700
ok.ids <- df[,max(winner), by=raceid][V1==1, raceid]
df <- df[raceid %in% ok.ids]
setkey(df, raceid, winodds)

## remove races for which no winner in the data eg: raceid 19700
ok.ids <- df[,max(winner), by=raceid][V1==1, raceid]
df <- df[raceid %in% ok.ids]
cara <- function(xs, theta) {
    (1 - exp(-theta*xs))/theta
}

## read 
df.maiden <- fread(paste(data.dir, "dataset_clm.csv", sep=""))
df.nonmaiden <- fread(paste(data.dir, "dataset_mcl.csv", sep=""))
df.salanie <- fread(paste(data.dir, "salanie_long.csv", sep=""))
setkey(df.maiden, raceid, winodds)
setkey(df.nonmaiden, raceid, winodds)
setkey(df.salanie, raceid, winodds)



solve.ps.st <- function(delta,theta, Rs, epsilon=10e-9) {
        n <- length(Rs)
        fun <- function(ps.plus.w) { # ps.plus.w = c(ps, w)
            ps <- ps.plus.w[1:n]
            w = ps.plus.w[n+1]
            d <- delta
            u <- cara 
            ## ws is d^(n-1) * ps * u(Rs, theta) + (d(1-d^n)/(1-d) - d^i)(1-ps)u(-1, theta)
            ds <- d^seq(1,n)
            rslt <- d^(n-1) * ps * u(Rs, theta) + (d*(1-d^n)/(1-d) - ds) * (1-ps) * u(-1, theta) - w
            return (c(rslt, sum(ps)-1))
        }
        ps <- rep(1/n, n) # initial value
        w <- 0.5
        rslt <- nleqslv(c(ps, w), fun)$x[1:n]
        rslt <- ifelse(rslt <= 0, epsilon, rslt)
        rslt
}


get.df.sample <- function(df, frac=0.1, random=FALSE, seed=10) {
    if (frac < 0 | frac > 1)
        stop("frac should be between 0 and 1")
    ids <- df[,unique(raceid)] 
    if(!random)
        set.seed(10)
    sids <- sample(ids, as.integer(frac*length(ids)))
    data.table(data.frame(df[raceid %in% sids]))
}

