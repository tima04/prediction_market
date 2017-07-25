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

## make clusters, cl is closed at the end of the file
NNODE = detectCores() 
cl <- makeForkCluster(NNODE)

df <- as.data.table(as.data.frame(df.salanie))

main <- function(N=3) {
    dt <- data.table(as.data.frame(df))
    train.size = floor(nrow(dt)/2)

    eu.score <- c()
    sal.score <- c()
    rdu.score <- c()

    for (i in 1:N) {
        set.seed(i)
        train.ind <- sample(1:nrow(dt), train.size)
        train.ind <- sort(train.ind)
        test.ind <- setdiff(1:nrow(dt), train.ind) 
        dt.train <- dt[train.ind]
        dt.test <- dt[test.ind]
        ##eu
        eu.par <- tryCatch(
            mle(minuslogl=function(t) -loglike.eu(t, dt.train),
                start=list(t=0.5) ,method="Brent", lower=-5, upper=5)@coef,
            error=function(e) {print("error in eu"); 0})
            

        eu.score <- c(eu.score, loglike.eu(eu.par, dt.test))
        
        ##sal
        sal.par <- tryCatch(
            mle(minuslogl=function(delta, t) -loglike.sal(delta, t, dt.train),
                       start=list(delta=0.5, t=-0.02),
                       method="Nelder-Mead")@coef, 
            error=function(e) {print("error in sal"); c(0, 0)})
        sal.score <- c(sal.score, loglike.sal(delta=sal.par[1], t=sal.par[2], dt.test))
        
        ##rdu
        rdu.par <- tryCatch(
            mle(minuslogl=function(a, beta, t) -loglike.rdu.parallel(a, beta, t, dt.train),
        ##rdu.par <- mle(minuslogl=function(a, beta, t) -loglike.rdu(a, beta, t, dt.train),
                start=list(a=0.5, beta=0.8, t=-0.02),
                method = "Nelder-Mead")@coef, 
            error=function(e) {print("error in sal"); c(0, 0, 0)})
        rdu.score <- c(rdu.score,
                       loglike.rdu.parallel(rdu.par[1], rdu.par[2], rdu.par[3], dt.test))
                       ## ifelse(all(rdu.par == c(0, 0, 0), NA,
                       ##            loglike.rdu.parallel(rdu.par[1], rdu.par[2], rdu.par[3], dt.test))))
    
        print(rdu.score)
    }


    cbind(eu.score, sal.score, rdu.score)
}

time <- system.time(rslt <- main(N=20))

try(stopCluster(cl), silent = TRUE)


foo <- function(x) {
    if (x < 0)
        stop("negative x")
    x^2
}

