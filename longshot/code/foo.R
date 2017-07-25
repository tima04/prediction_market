pi.l <- function(x, l, theta=0.5) {
    (2*x + theta)/(2*(2*x + theta - l))
}


pi.g <- function(x, l, theta=0.5) {
    1 - pi.l(x, l, theta)
}

gain <- function(l, pi.l) {
    l*pi.l/(1-pi.l)
}

x <- 100
l <- 80
pi.star <- pi.g(x, l, 0.5) 
g <- gain(l, 1-pi.star)
pi.star
g

wght <- function(pi, delta=0.5) {
    pi/(pi + (1-pi)*delta)
}

premium <- function(x, l, pi.l, delta=0.5) {
    pi.g <- 1 - pi.l
    pig1 <- pi.g/(pi.g + (1-pi.g)*delta) 
    v1 <- pig1*(x + gain(l, pi.l)) + (1-pig1)*(x - l)
    x - v1
}

wgt.prelec <- function(p, alpha=0.63) {
    exp(-(-log(p))^alpha)
}

premium.pt <- function(x, l, pi.l, alpha=0.5) {
    pi.g <- 1-pi.l
    wgt <- wgt.prelec(pi.g)
    pt <- wgt * (x + gain(l, pi.l)) + (1 - wgt) * (x - l) 
    x - pt
}

x <- 100
l <- 10
pi.l <- .99999
gain(l, pi.l)
premium(x, l, pi.l)
premium.pt(x, l, pi.l)

delta <- 0.1 
ps <- seq(0,0.2,0.00001)
ws.pt <- wgt.prelec(ps, 0.8)/ps  
ws.sl <- ps/(ps + (1-ps)*delta)/ps 
plot(ws.pt~ps, type="l",log="y")
lines(ws.sl~ps, col="red")
