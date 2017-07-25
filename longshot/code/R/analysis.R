library(data.table)
data.dir <- "../data/gandhi_restd/" 

df.clm <- data.table(read.csv(paste(data.dir, "dataset_clm.csv", sep=""))) ## maiden horses
df.mcl <- data.table(read.csv(paste(data.dir, "dataset_mcl.csv", sep=""))) ## non-maiden(experienced) 

## combining 
df.clm[, maiden:=1]
df.mcl[, maiden:=0]
df <- as.data.table(rbind(df.clm, df.mcl))
setkey(df, raceid, winodds)

## calculating win probablility associated with each horse, from bet data
df[ ,total.raw.prob := sum(1/(1+winodds)), by=raceid]
df[, prob := 1/(1+winodds)/total.raw.prob]

## win probablitites of favorite

##======================================== playground

num.races <- df[, length(unique(raceid))] # =  116397
n = 7  
xs <- c()
odds <- c()
ps <- c()
is.maiden <- 0 
max.win.odd <- 50   
races <- df[numhorses == n & (maiden==is.maiden), max(winodds), by=raceid][V1 > max.win.odd, raceid]
for (i in 1:n) {
    xs <- c(xs,
            df[(numhorses == n) & (maiden==is.maiden) & (raceid %in% races), winodds[i]*winner[i], by=raceid][, sum(V1)/.N])
    odds <- c(odds,
            df[(numhorses == n) & (maiden==is.maiden) & (raceid %in% races), winodds[i], by=raceid][, median(V1)])
    ps <- c(ps,
            df[(numhorses == n) & (maiden==is.maiden) & (raceid %in% races), winner[i], by=raceid][, sum(V1)/.N])
    }
print("============")
xs
print(paste(which.max(xs), n, sep=", "))
ps
odds

## first, second, second-last, last
xs <- matrix(ncol=11, nrow=4) 
odds <- matrix(ncol=11, nrow=4) 
ps <- matrix(ncol=11, nrow=4) 
df2 <- df[numhorses > 3]
is.maiden <- 0 
max.win.odd <- 40   
#races <- df2[(maiden==is.maiden), max(winodds), by=raceid][V1 > max.win.odd, raceid]
races <- df2[, max(winodds), by=raceid][V1 > max.win.odd, raceid]
for (n in 4:14) {
    for (i in 1:4) {
        j <- ifelse(i < 3, i, ifelse(i==3, n-1, n))
        xs[i,n-3] <- df2[(numhorses == n) & (maiden==is.maiden) & (raceid %in% races), (1+winodds[j])*winner[j], by=raceid][, sum(V1)/.N]
        odds[i,n-3] <- df2[(numhorses == n) & (maiden==is.maiden) & (raceid %in% races), winodds[j], by=raceid][, median(V1)]
        ps[i,n-3] <- df2[(numhorses == n) & (maiden==is.maiden) & (raceid %in% races), winner[j], by=raceid][, sum(V1)/.N]
    }
}
print("============")
#apply(xs, 1, median, na.rm=T)
apply(ps, 1, median, na.rm=T)
apply(odds, 1, median, na.rm=T)
xs[,3:7 ]

df[maiden == 0 & winner==1 & winodds > 40, .N]
ids <- df[maiden == 0 & winner==1 & winodds > 40, raceid] 
df[raceid == ids[i]-1, winodds]
df[raceid == ids[i], winodds]
df[raceid == ids[i]+1, winodds]

rslt <- matrix(nrow=10, ncol=2)
for(i in 1:10) {
    ali <- df[maiden==0 & numhorses==10]
    ali[, length(unique(raceid))]
    foo <- ali[, .(sub=prob[i], obj=winner[i]), by=raceid][, .(sub=mean(sub), obj=mean(obj))]
    rslt[i, 1] <- foo[,sub]
    rslt[i, 2] <- foo[,obj]
}

ali.rslt <- data.table(rslt)
names(ali.rslt) <- c("sub", "obj")
N <- ali[,.N]
ali.rslt[, se:=sqrt((obj*(1-obj)/N))]
ali.rslt[, cv := 100*se/obj]
ali.rslt[, diff := (sub-obj)/se]
    
