library(data.table)
data.dir <- "../data/gandhi_restd/" 

df.clm <- data.table(read.csv(paste(data.dir, "dataset_clm.csv", sep=""))) ## maiden horses
df.mcl <- data.table(read.csv(paste(data.dir, "dataset_mcl.csv", sep=""))) ## non-maiden(experienced) 

## combining 
df.clm[, maiden:=1]
df.mcl[, maiden:=0]
df <- as.data.table(rbind(df.clm, df.mcl))
setkey(df, raceid, winodds)

## get the table 2 of Ali from the gandhi data.
## args: mdn is 0 or 1, 1 for maiden horses and 0  for non-maiden.
get.ali.rslt <- function(mdn=0) {
    ## calculating win probablility associated with each horse, from bet data
    df[ ,total.raw.prob := sum(1/(1+winodds)), by=raceid]
    df[, prob := 1/(1+winodds)/total.raw.prob]

    rslt <- matrix(nrow=10, ncol=4)
    ali <- df[maiden==mdn & numhorses %in% 6:10]
    for(i in 1:10) {
        foo <- ali[, .(sub=prob[i], obj=winner[i], odd=winodds[i]), by=raceid][, .(sub=mean(sub, na.rm=T),
                                                                                   obj=mean(obj, na.rm=T),
                                                                                   odd=mean(odd, na.rm=T))]
        rslt[i, 1] <- foo[,sub]
        rslt[i, 2] <- foo[,obj]
        rslt[i, 3] <- ali[numhorses >= i, length(unique(raceid))] 
        rslt[i, 4] <- foo[,odd]
    }
    ali.rslt <- data.table(rslt)
    names(ali.rslt) <- c("sub", "obj", "N", "odd")

    ali.rslt[, se:=sqrt(obj*(1-obj)/N)]
    ali.rslt[, ":=" (cv=100*se/obj, diff=(sub-obj)/se)]
    ali.rslt <- data.table(apply(ali.rslt, 2, round, digits=3))
    ali.rslt
}

## regression 
ali.rslt[, ":=" (x = (1 + odd), u = obj[8]/obj)]
fit <- lm(log10(u) ~ log10(x), data=ali.rslt[1:8])
summary(fit)

## table2 ali (page 808)
ali.table <- data.table(
    sub=c(.3237, .2077, .1513, .1121, .0827, .0601, .0417, .0276, NA, NA),
    obj=c(.3583, .2049, .1526, .1047, .0762, .0552, .0341, .0206, .0033, .0141),
    N=c(rep(20247, 4), 20231, 20088, 19281, 15749, 299, 71))
ali.table[, se:=sqrt(obj*(1-obj)/N)]
ali.table[,":=" (cv=100*se/obj, diff=(sub-obj)/se)]
## table 3
ali.table[, ":=" (x = c(2.5, 4, 5.4, 7.2, 9.7, 13.3, 18.7, 27.4, NA, NA))]
ali.table[, ":=" (u = obj[8]/obj)]
fit <- lm(log10(u) ~ log10(x), data=ali.table[1:8])
summary(fit)
plot(u~x, data=ali.table[1:8], type="b")
