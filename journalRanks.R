## Journal Ranks analysis
## Bradshaw & Brook
## Sep 2014/updated Mar 2015/updated Feb 2016
## Added to GitHub 29/06/2017

## Remove everything
rm(list = ls())

## libraries
library(Hmisc)
library(cluster)
library(bootstrap)

## set working directory
setwd("~/Documents/Papers/Other/Journal Ranks/data/")

###############################################
## import journal lists
dat.ecol <- read.table("ecol.csv", header=T, sep=",")
dat.ecol.samp <- read.table("ecolsamp.csv", header=T, sep=",")
dat.ecol.samp14 <- read.table("ecolsamp2014.csv", header=T, sep=",")
dat.ecol.samp15 <- read.table("ecolsamp2015.csv", header=T, sep=",")
dat.medi <- read.table("medi.csv", header=T, sep=",")
dat.mult <- read.table("mult.csv", header=T, sep=",")
dat.mish <- read.table("marfish.csv", header=T, sep=",")
dat.ogyn <- read.table("obsgyn.csv", header=T, sep=",") 
dat.ener <- read.table("energy.csv", header=T, sep=",") 
lecol <- dim(dat.ecol)[1]
lmedi <- dim(dat.medi)[1]
lmult <- dim(dat.mult)[1]
lmish <- dim(dat.mish)[1]
logyn <- dim(dat.ogyn)[1]
lener <- dim(dat.ener)[1]

## 'BIOLOGY' journals listed in Web of Science (correlation between IF, EFS, AIS)
dat.biol <- read.table("biology.csv", header=T, sep=",")
biol.cor <- dat.biol[,c(4,8,9)]
biol.cormat <- cor(biol.cor,y=NULL,use="complete.obs",method="spearman")
lvar <- dim(biol.cormat)[2]
for (c in 1:lvar) {
  biol.cormat[c,c:lvar] <- NA}
round(biol.cormat[,-lvar], 3)


## Ecology
## Plots
par(mfrow=c(3,2))
plot((dat.ecol$n), dat.ecol$h5, pch=19, xlab="articles", ylab="h5")
plot((dat.ecol$h5), dat.ecol$h5med, pch=19, xlab="h5", ylab="h5med")
plot((dat.ecol$h5), dat.ecol$IF2013, pch=19, xlab="h5", ylab="IF2013")
plot((dat.ecol$IF2013), dat.ecol$IF5, pch=19, xlab="IF2013", ylab="IF5")
plot((dat.ecol$IF2013), dat.ecol$IM, pch=19, xlab="IF2013", ylab="Immediacy")
plot((dat.ecol$h5/dat.ecol$n), dat.ecol$IF2013, pch=19, xlab="h5/article", ylab="IF2013")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot((dat.ecol$SNIP13), dat.ecol$IPP13, pch=19, xlab="SNIP13", ylab="IPP13")
plot((dat.ecol$IPP13), dat.ecol$SJR13, pch=19, xlab="IPP13", ylab="SJR13")
plot((dat.ecol$SNIP13), dat.ecol$SJR13, pch=19, xlab="SNIP13", ylab="SJR13")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot((dat.ecol$IF5), dat.ecol$IPP13, pch=19, xlab="IF5", ylab="IPP13")
plot((dat.ecol$h5), dat.ecol$SJR13, pch=19, xlab="h5", ylab="SJR13")
plot((dat.ecol$h5med), dat.ecol$SJR13, pch=19, xlab="h5med", ylab="SJR13")
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot((dat.ecol$n), dat.ecol$IF2013, pch=19, xlab="articles", ylab="IF2013")
plot((dat.ecol$cites/dat.ecol$n), dat.ecol$IM, pch=19, xlab="cites/article", ylab="IM")
par(mfrow=c(1,1))

## Ranks
cites.rnk <- lecol - rank(dat.ecol$cites, ties.method="average") + 1
h5n.rnk <- lecol - rank(dat.ecol$h5/dat.ecol$n, ties.method="average") + 1
h5med.rnk <- lecol - rank(dat.ecol$h5med, ties.method="average") + 1
IF2013.rnk <- lecol - rank(dat.ecol$IF2013, ties.method="average") + 1
IM.rnk <- lecol - rank(dat.ecol$IM, ties.method="average") + 1
SNIP.rnk <- lecol - rank(dat.ecol$SNIP13, ties.method="average") + 1
IPP.rnk <- lecol - rank(dat.ecol$IPP13, ties.method="average") + 1
SJR.rnk <- lecol - rank(dat.ecol$SJR13, ties.method="average") + 1

plot(cites.rnk, IF2013.rnk, pch=19, xlab="citations", ylab="IF5")
fit <- lm(IF2013.rnk ~ cites.rnk)
abline(fit)

## data frame
ecol.rnk.dat <- data.frame(dat.ecol$Journal, h5n.rnk, h5med.rnk, IF2013.rnk, IM.rnk, SNIP.rnk, IPP.rnk, SJR.rnk)
colnames(ecol.rnk.dat) <- c("journal","h5n", "h5med", "IF2013", "IM", "SNIP", "IPP", "SJR")

## Average Rank
ecol.avg.rnk <- apply(ecol.rnk.dat[,2:8], 1, "mean")
ecol.sd.rnk <- apply(ecol.rnk.dat[,2:8], 1, "sd")
ecol.cv.rnk <- 100*ecol.sd.rnk/ecol.avg.rnk
ecol.rnk.dat <- data.frame(ecol.rnk.dat, ecol.avg.rnk, ecol.sd.rnk, ecol.cv.rnk)

ecol.rnk.sort <- ecol.rnk.dat[order(ecol.rnk.dat[,9],decreasing=F),]
ecol.rnk.sort

plot(ecol.rnk.sort$ecol.avg.rnk,ecol.rnk.sort$ecol.cv.rnk, pch=19)
plot(ecol.rnk.sort$ecol.avg.rnk,ecol.rnk.sort$ecol.sd.rnk, pch=19)

# Histograms
hist(ecol.rnk.sort$ecol.avg.rnk)


## Medical
## Plots
par(mfrow=c(3,2))
plot((dat.medi$n), dat.medi$h5, pch=19, xlab="articles", ylab="h5")
plot((dat.medi$h5), dat.medi$h5med, pch=19, xlab="h5", ylab="h5med")
plot((dat.medi$h5), dat.medi$IF2013, pch=19, xlab="h5", ylab="IF2013")
plot((dat.medi$IF2013), dat.medi$IF5, pch=19, xlab="IF2013", ylab="IF5")
plot((dat.medi$IF2013), dat.medi$IM, pch=19, xlab="IF2013", ylab="Immediacy")
plot((dat.medi$h5/dat.medi$n), dat.medi$IF2013, pch=19, xlab="h5/article", ylab="IF2013")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot((dat.medi$SNIP13), dat.medi$IPP13, pch=19, xlab="SNIP13", ylab="IPP13")
plot((dat.medi$IPP13), dat.medi$SJR13, pch=19, xlab="IPP13", ylab="SJR13")
plot((dat.medi$SNIP13), dat.medi$SJR13, pch=19, xlab="SNIP13", ylab="SJR13")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot((dat.medi$IF5), dat.medi$IPP13, pch=19, xlab="IF5", ylab="IPP13")
plot((dat.medi$h5), dat.medi$SJR13, pch=19, xlab="h5", ylab="SJR13")
plot((dat.medi$h5med), dat.medi$SJR13, pch=19, xlab="h5med", ylab="SJR13")
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot((dat.medi$n), dat.medi$IF2013, pch=19, xlab="articles", ylab="IF2013")
plot((dat.medi$cites/dat.medi$n), dat.medi$IM, pch=19, xlab="cites/article", ylab="IM")
par(mfrow=c(1,1))

## Ranks
cites.rnk <- lmedi - rank(dat.medi$cites, ties.method="average") + 1
h5n.rnk <- lmedi - rank(dat.medi$h5/dat.medi$n, ties.method="average") + 1
h5med.rnk <- lmedi - rank(dat.medi$h5med, ties.method="average") + 1
IF2013.rnk <- lmedi - rank(dat.medi$IF2013, ties.method="average") + 1
IM.rnk <- lmedi - rank(dat.medi$IM, ties.method="average") + 1
SNIP.rnk <- lmedi - rank(dat.medi$SNIP13, ties.method="average") + 1
IPP.rnk <- lmedi - rank(dat.medi$IPP13, ties.method="average") + 1
SJR.rnk <- lmedi - rank(dat.medi$SJR13, ties.method="average") + 1

plot(cites.rnk, IF2013.rnk, pch=19, xlab="citations", ylab="IF5")
fit <- lm(IF2013.rnk ~ cites.rnk)
abline(fit)


## data frame
medi.rnk.dat <- data.frame(dat.medi$Journal, h5n.rnk, h5med.rnk, IF2013.rnk, IM.rnk, SNIP.rnk, IPP.rnk, SJR.rnk)
colnames(medi.rnk.dat) <- c("journal","h5n", "h5med", "IF2013", "IM", "SNIP", "IPP", "SJR")

## Average Rank
medi.avg.rnk <- apply(medi.rnk.dat[,2:8], 1, "mean")
medi.sd.rnk <- apply(medi.rnk.dat[,2:8], 1, "sd")
medi.cv.rnk <- 100*medi.sd.rnk/medi.avg.rnk
medi.rnk.dat <- data.frame(medi.rnk.dat, medi.avg.rnk, medi.sd.rnk, medi.cv.rnk)

medi.rnk.sort <- medi.rnk.dat[order(medi.rnk.dat[,9],decreasing=F),]
medi.rnk.sort

plot(medi.rnk.sort$medi.avg.rnk,medi.rnk.sort$medi.cv.rnk, pch=19)
plot(medi.rnk.sort$medi.avg.rnk,medi.rnk.sort$medi.sd.rnk, pch=19)

# Histograms
hist(medi.rnk.sort$medi.avg.rnk)


## Multidisciplinary
## Plots
par(mfrow=c(3,2))
plot((dat.mult$n), dat.mult$h5, pch=19, xlab="articles", ylab="h5")
plot((dat.mult$h5), dat.mult$h5med, pch=19, xlab="h5", ylab="h5med")
plot((dat.mult$h5), dat.mult$IF2013, pch=19, xlab="h5", ylab="IF2013")
plot((dat.mult$IF2013), dat.mult$IF5, pch=19, xlab="IF2013", ylab="IF5")
plot((dat.mult$IF2013), dat.mult$IM, pch=19, xlab="IF2013", ylab="Immultacy")
plot((dat.mult$h5/dat.mult$n), dat.mult$IF2013, pch=19, xlab="h5/article", ylab="IF2013")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot((dat.mult$SNIP13), dat.mult$IPP13, pch=19, xlab="SNIP13", ylab="IPP13")
plot((dat.mult$IPP13), dat.mult$SJR13, pch=19, xlab="IPP13", ylab="SJR13")
plot((dat.mult$SNIP13), dat.mult$SJR13, pch=19, xlab="SNIP13", ylab="SJR13")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot((dat.mult$IF5), dat.mult$IPP13, pch=19, xlab="IF5", ylab="IPP13")
plot((dat.mult$h5), dat.mult$SJR13, pch=19, xlab="h5", ylab="SJR13")
plot((dat.mult$h5med), dat.mult$SJR13, pch=19, xlab="h5med", ylab="SJR13")
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot((dat.mult$n), dat.mult$IF2013, pch=19, xlab="articles", ylab="IF2013")
plot((dat.mult$cites/dat.mult$n), dat.mult$IM, pch=19, xlab="cites/article", ylab="IM")
par(mfrow=c(1,1))

## Ranks
cites.rnk <- lmult - rank(dat.mult$cites, ties.method="average") + 1
h5n.rnk <- lmult - rank(dat.mult$h5/dat.mult$n, ties.method="average") + 1
h5med.rnk <- lmult - rank(dat.mult$h5med, ties.method="average") + 1
IF2013.rnk <- lmult - rank(dat.mult$IF2013, ties.method="average") + 1
IM.rnk <- lmult - rank(dat.mult$IM, ties.method="average") + 1
SNIP.rnk <- lmult - rank(dat.mult$SNIP13, ties.method="average") + 1
IPP.rnk <- lmult - rank(dat.mult$IPP13, ties.method="average") + 1
SJR.rnk <- lmult - rank(dat.mult$SJR13, ties.method="average") + 1

plot(cites.rnk, IF2013.rnk, pch=19, xlab="citations", ylab="IF5")
fit <- lm(IF2013.rnk ~ cites.rnk)
abline(fit)

## data frame
mult.rnk.dat <- data.frame(dat.mult$Journal, h5n.rnk, h5med.rnk, IF2013.rnk, IM.rnk, SNIP.rnk, IPP.rnk, SJR.rnk)
colnames(mult.rnk.dat) <- c("journal","h5n", "h5med", "IF2013", "IM", "SNIP", "IPP", "SJR")

## Average Rank
mult.avg.rnk <- apply(mult.rnk.dat[,2:8], 1, "mean")
mult.sd.rnk <- apply(mult.rnk.dat[,2:8], 1, "sd")
mult.cv.rnk <- 100*mult.sd.rnk/mult.avg.rnk
mult.rnk.dat <- data.frame(mult.rnk.dat, mult.avg.rnk, mult.sd.rnk, mult.cv.rnk)

mult.rnk.sort <- mult.rnk.dat[order(mult.rnk.dat[,9],decreasing=F),]
mult.rnk.sort

plot(mult.rnk.sort$mult.avg.rnk,mult.rnk.sort$mult.cv.rnk, pch=19)
plot(mult.rnk.sort$mult.avg.rnk,mult.rnk.sort$mult.sd.rnk, pch=19)

# Histograms
hist(mult.rnk.sort$mult.avg.rnk)



## Marine & Fisheries
## Plots
par(mfrow=c(3,2))
plot((dat.mish$n), dat.mish$h5, pch=19, xlab="articles", ylab="h5")
plot((dat.mish$h5), dat.mish$h5med, pch=19, xlab="h5", ylab="h5med")
plot((dat.mish$h5), dat.mish$IF2013, pch=19, xlab="h5", ylab="IF2013")
plot((dat.mish$IF2013), dat.mish$IF5, pch=19, xlab="IF2013", ylab="IF5")
plot((dat.mish$IF2013), dat.mish$IM, pch=19, xlab="IF2013", ylab="Immishacy")
plot((dat.mish$h5/dat.mish$n), dat.mish$IF2013, pch=19, xlab="h5/article", ylab="IF2013")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot((dat.mish$SNIP13), dat.mish$IPP13, pch=19, xlab="SNIP13", ylab="IPP13")
plot((dat.mish$IPP13), dat.mish$SJR13, pch=19, xlab="IPP13", ylab="SJR13")
plot((dat.mish$SNIP13), dat.mish$SJR13, pch=19, xlab="SNIP13", ylab="SJR13")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot((dat.mish$IF5), dat.mish$IPP13, pch=19, xlab="IF5", ylab="IPP13")
plot((dat.mish$h5), dat.mish$SJR13, pch=19, xlab="h5", ylab="SJR13")
plot((dat.mish$h5med), dat.mish$SJR13, pch=19, xlab="h5med", ylab="SJR13")
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot((dat.mish$n), dat.mish$IF2013, pch=19, xlab="articles", ylab="IF2013")
plot((dat.mish$cites/dat.mish$n), dat.mish$IM, pch=19, xlab="cites/article", ylab="IM")
par(mfrow=c(1,1))

## Ranks
cites.rnk <- lmish - rank(dat.mish$cites, ties.method="average") + 1
h5n.rnk <- lmish - rank(dat.mish$h5/dat.mish$n, ties.method="average") + 1
h5med.rnk <- lmish - rank(dat.mish$h5med, ties.method="average") + 1
IF2013.rnk <- lmish - rank(dat.mish$IF2013, ties.method="average") + 1
IM.rnk <- lmish - rank(dat.mish$IM, ties.method="average") + 1
SNIP.rnk <- lmish - rank(dat.mish$SNIP13, ties.method="average") + 1
IPP.rnk <- lmish - rank(dat.mish$IPP13, ties.method="average") + 1
SJR.rnk <- lmish - rank(dat.mish$SJR13, ties.method="average") + 1

plot(cites.rnk, IF2013.rnk, pch=19, xlab="citations", ylab="IF5")
fit <- lm(IF2013.rnk ~ cites.rnk)
abline(fit)

## data frame
mish.rnk.dat <- data.frame(dat.mish$Journal, h5n.rnk, h5med.rnk, IF2013.rnk, IM.rnk, SNIP.rnk, IPP.rnk, SJR.rnk)
colnames(mish.rnk.dat) <- c("journal","h5n", "h5med", "IF2013", "IM", "SNIP", "IPP", "SJR")

## Average Rank
mish.avg.rnk <- apply(mish.rnk.dat[,2:8], 1, "mean")
mish.sd.rnk <- apply(mish.rnk.dat[,2:8], 1, "sd")
mish.cv.rnk <- 100*mish.sd.rnk/mish.avg.rnk
mish.rnk.dat <- data.frame(mish.rnk.dat, mish.avg.rnk, mish.sd.rnk, mish.cv.rnk)

mish.rnk.sort <- mish.rnk.dat[order(mish.rnk.dat[,9],decreasing=F),]
mish.rnk.sort

plot(mish.rnk.sort$mish.avg.rnk,mish.rnk.sort$mish.cv.rnk, pch=19)
plot(mish.rnk.sort$mish.avg.rnk,mish.rnk.sort$mish.sd.rnk, pch=19)

# Histograms
hist(mish.rnk.sort$mish.avg.rnk)


## Obstetrics & Gynecology
## Plots
par(mfrow=c(3,2))
plot((dat.ogyn$n), dat.ogyn$h5, pch=19, xlab="articles", ylab="h5")
plot((dat.ogyn$h5), dat.ogyn$h5med, pch=19, xlab="h5", ylab="h5med")
plot((dat.ogyn$h5), dat.ogyn$IF2013, pch=19, xlab="h5", ylab="IF2013")
plot((dat.ogyn$IF2013), dat.ogyn$IF5, pch=19, xlab="IF2013", ylab="IF5")
plot((dat.ogyn$IF2013), dat.ogyn$IM, pch=19, xlab="IF2013", ylab="Imogynacy")
plot((dat.ogyn$h5/dat.ogyn$n), dat.ogyn$IF2013, pch=19, xlab="h5/article", ylab="IF2013")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot((dat.ogyn$SNIP13), dat.ogyn$IPP13, pch=19, xlab="SNIP13", ylab="IPP13")
plot((dat.ogyn$IPP13), dat.ogyn$SJR13, pch=19, xlab="IPP13", ylab="SJR13")
plot((dat.ogyn$SNIP13), dat.ogyn$SJR13, pch=19, xlab="SNIP13", ylab="SJR13")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
plot((dat.ogyn$IF5), dat.ogyn$IPP13, pch=19, xlab="IF5", ylab="IPP13")
plot((dat.ogyn$h5), dat.ogyn$SJR13, pch=19, xlab="h5", ylab="SJR13")
plot((dat.ogyn$h5med), dat.ogyn$SJR13, pch=19, xlab="h5med", ylab="SJR13")
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot((dat.ogyn$n), dat.ogyn$IF2013, pch=19, xlab="articles", ylab="IF2013")
plot((dat.ogyn$cites/dat.ogyn$n), dat.ogyn$IM, pch=19, xlab="cites/article", ylab="IM")
par(mfrow=c(1,1))

## Ranks
cites.rnk <- logyn - rank(dat.ogyn$cites, ties.method="average") + 1
h5n.rnk <- logyn - rank(dat.ogyn$h5/dat.ogyn$n, ties.method="average") + 1
h5med.rnk <- logyn - rank(dat.ogyn$h5med, ties.method="average") + 1
IF2013.rnk <- logyn - rank(dat.ogyn$IF2013, ties.method="average") + 1
IM.rnk <- logyn - rank(dat.ogyn$IM, ties.method="average") + 1
SNIP.rnk <- logyn - rank(dat.ogyn$SNIP13, ties.method="average") + 1
IPP.rnk <- logyn - rank(dat.ogyn$IPP13, ties.method="average") + 1
SJR.rnk <- logyn - rank(dat.ogyn$SJR13, ties.method="average") + 1

plot(cites.rnk, IF2013.rnk, pch=19, xlab="citations", ylab="IF5")
fit <- lm(IF2013.rnk ~ cites.rnk)
abline(fit)

## data frame
ogyn.rnk.dat <- data.frame(dat.ogyn$Journal, h5n.rnk, h5med.rnk, IF2013.rnk, IM.rnk, SNIP.rnk, IPP.rnk, SJR.rnk)
colnames(ogyn.rnk.dat) <- c("journal","h5n", "h5med", "IF2013", "IM", "SNIP", "IPP", "SJR")

## Average Rank
ogyn.avg.rnk <- apply(ogyn.rnk.dat[,2:8], 1, "mean")
ogyn.sd.rnk <- apply(ogyn.rnk.dat[,2:8], 1, "sd")
ogyn.cv.rnk <- 100*ogyn.sd.rnk/ogyn.avg.rnk
ogyn.rnk.dat <- data.frame(ogyn.rnk.dat, ogyn.avg.rnk, ogyn.sd.rnk, ogyn.cv.rnk)

ogyn.rnk.sort <- ogyn.rnk.dat[order(ogyn.rnk.dat[,9],decreasing=F),]
ogyn.rnk.sort

plot(ogyn.rnk.sort$ogyn.avg.rnk,ogyn.rnk.sort$ogyn.cv.rnk, pch=19)
plot(ogyn.rnk.sort$ogyn.avg.rnk,ogyn.rnk.sort$ogyn.sd.rnk, pch=19)

# Histograms
hist(ogyn.rnk.sort$ogyn.avg.rnk)



## Correlation matrices
## correlation matrix for GLM input data
ecol.cor <- dat.ecol[,3:12]
ecol.cormat <- cor(ecol.cor,y=NULL,use="complete.obs",method="spearman")
lvar <- dim(ecol.cormat)[2]
for (c in 1:lvar) {
  ecol.cormat[c,c:lvar] <- NA}
round(ecol.cormat[,-lvar], 3)

medi.cor <- dat.medi[,3:12]
medi.cormat <- cor(medi.cor,y=NULL,use="complete.obs",method="spearman")
lvar <- dim(medi.cormat)[2]
for (c in 1:lvar) {
  medi.cormat[c,c:lvar] <- NA}
round(medi.cormat[,-lvar], 3)

mult.cor <- dat.mult[,3:12]
mult.cormat <- cor(mult.cor,y=NULL,use="complete.obs",method="spearman")
lvar <- dim(mult.cormat)[2]
for (c in 1:lvar) {
  mult.cormat[c,c:lvar] <- NA}
round(mult.cormat[,-lvar], 3)

mish.cor <- dat.mish[,3:12]
mish.cormat <- cor(mish.cor,y=NULL,use="complete.obs",method="spearman")
lvar <- dim(mish.cormat)[2]
for (c in 1:lvar) {
  mish.cormat[c,c:lvar] <- NA}
round(mish.cormat[,-lvar], 3)

ogyn.cor <- dat.ogyn[,3:12]
ogyn.cormat <- cor(ogyn.cor,y=NULL,use="complete.obs",method="spearman")
lvar <- dim(ogyn.cormat)[2]
for (c in 1:lvar) {
  ogyn.cormat[c,c:lvar] <- NA}
round(ogyn.cormat[,-lvar], 3)


## Combined plots & analysis
par(mfrow=c(1,3))
plot(log10(dat.ecol$n), dat.ecol$h5, pch=19, xlab="log n articles", ylab="h5")
plot(log10(dat.medi$n), dat.medi$h5, pch=19, xlab="log n articles", ylab="h5")
plot(log10(dat.mult$n), dat.mult$h5, pch=19, xlab="log n articles", ylab="h5")
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot(log10(dat.mish$n), dat.mish$h5, pch=19, xlab="log n articles", ylab="h5")
plot(log10(dat.ogyn$n), dat.ogyn$h5, pch=19, xlab="log n articles", ylab="h5")
par(mfrow=c(1,1))


attach(dat.ecol)
h5ln.ecol <- h5/log10(n)
detach(dat.ecol)

attach(dat.medi)
h5ln.medi <- h5/log10(n)
detach(dat.medi)

attach(dat.mult)
h5ln.mult <- h5/log10(n)
detach(dat.mult)

attach(dat.mish)
h5ln.mish <- h5/log10(n)
detach(dat.mish)

attach(dat.ogyn)
h5ln.ogyn <- h5/log10(n)
detach(dat.ogyn)

attach(dat.ener)
h5ln.ener <- h5/log10(n)
detach(dat.ener)

par(mfrow=c(1,3))
plot(h5ln.ecol, dat.ecol$IF2013, pch=19, xlab="h5/log n articles", ylab="IF2013")
plot(h5ln.medi, dat.medi$IF2013, pch=19, xlab="h5/log n articles", ylab="IF2013")
plot(h5ln.mult, dat.mult$IF2013, pch=19, xlab="h5/log n articles", ylab="IF2013")
par(mfrow=c(1,1))


# final set
ecol.h5ln.rnk <- lecol - rank(h5ln.ecol, ties.method="average") + 1
ecol.IF13.rnk <- lecol - rank(dat.ecol$IF2013, ties.method="average") + 1
ecol.IM.rnk <- lecol - rank(dat.ecol$IM, ties.method="average") + 1
ecol.SNIP13.rnk <- lecol - rank(dat.ecol$SNIP13, ties.method="average") + 1
ecol.SJR13.rnk <- lecol - rank(dat.ecol$SJR13, ties.method="average") + 1
ecol.rnk.dat <- data.frame(dat.ecol$Journal, ecol.h5ln.rnk, ecol.IF13.rnk, ecol.IM.rnk, ecol.SNIP13.rnk, ecol.SJR13.rnk)
colnames(ecol.rnk.dat) <- c("journal","h5ln", "IF13", "IM", "SNIP", "SJR")
ecol.avg.rnk <- apply(ecol.rnk.dat[,2:6], 1, "mean")
ecol.sd.rnk <- apply(ecol.rnk.dat[,2:6], 1, "sd")
ecol.rnk.dat <- data.frame(ecol.rnk.dat, ecol.avg.rnk, ecol.sd.rnk)
ecol.rnk.sort <- ecol.rnk.dat[order(ecol.rnk.dat[,7],decreasing=F),]
ecol.rnk.sort
plot(ecol.rnk.sort$ecol.avg.rnk, ecol.rnk.sort$ecol.sd.rnk)

medi.h5ln.rnk <- lmedi - rank(h5ln.medi, ties.method="average") + 1
medi.IF13.rnk <- lmedi - rank(dat.medi$IF2013, ties.method="average") + 1
medi.IM.rnk <- lmedi - rank(dat.medi$IM, ties.method="average") + 1
medi.SNIP13.rnk <- lmedi - rank(dat.medi$SNIP13, ties.method="average") + 1
medi.SJR13.rnk <- lmedi - rank(dat.medi$SJR13, ties.method="average") + 1
medi.rnk.dat <- data.frame(dat.medi$Journal, medi.h5ln.rnk, medi.IF13.rnk, medi.IM.rnk, medi.SNIP13.rnk, medi.SJR13.rnk)
colnames(medi.rnk.dat) <- c("journal","h5ln", "IF13", "IM", "SNIP", "SJR")
medi.avg.rnk <- apply(medi.rnk.dat[,2:6], 1, "mean")
medi.sd.rnk <- apply(medi.rnk.dat[,2:6], 1, "sd")
medi.rnk.dat <- data.frame(medi.rnk.dat, medi.avg.rnk, medi.sd.rnk)
medi.rnk.sort <- medi.rnk.dat[order(medi.rnk.dat[,7],decreasing=F),]
medi.rnk.sort
plot(medi.rnk.sort$medi.avg.rnk, medi.rnk.sort$medi.sd.rnk)

mult.h5ln.rnk <- lmult - rank(h5ln.mult, ties.method="average") + 1
mult.IF13.rnk <- lmult - rank(dat.mult$IF2013, ties.method="average") + 1
mult.IM.rnk <- lmult - rank(dat.mult$IM, ties.method="average") + 1
mult.SNIP13.rnk <- lmult - rank(dat.mult$SNIP13, ties.method="average") + 1
mult.SJR13.rnk <- lmult - rank(dat.mult$SJR13, ties.method="average") + 1
mult.rnk.dat <- data.frame(dat.mult$Journal, mult.h5ln.rnk, mult.IF13.rnk, mult.IM.rnk, mult.SNIP13.rnk, mult.SJR13.rnk)
colnames(mult.rnk.dat) <- c("journal","h5ln", "IF13", "IM", "SNIP", "SJR")
mult.avg.rnk <- apply(mult.rnk.dat[,2:6], 1, "mean")
mult.sd.rnk <- apply(mult.rnk.dat[,2:6], 1, "sd")
mult.rnk.dat <- data.frame(mult.rnk.dat, mult.avg.rnk, mult.sd.rnk)
mult.rnk.sort <- mult.rnk.dat[order(mult.rnk.dat[,7],decreasing=F),]
mult.rnk.sort
plot(mult.rnk.sort$mult.avg.rnk, mult.rnk.sort$mult.sd.rnk)

mish.h5ln.rnk <- lmish - rank(h5ln.mish, ties.method="average") + 1
mish.IF13.rnk <- lmish - rank(dat.mish$IF2013, ties.method="average") + 1
mish.IM.rnk <- lmish - rank(dat.mish$IM, ties.method="average") + 1
mish.SNIP13.rnk <- lmish - rank(dat.mish$SNIP13, ties.method="average") + 1
mish.SJR13.rnk <- lmish - rank(dat.mish$SJR13, ties.method="average") + 1
mish.rnk.dat <- data.frame(dat.mish$Journal, mish.h5ln.rnk, mish.IF13.rnk, mish.IM.rnk, mish.SNIP13.rnk, mish.SJR13.rnk)
colnames(mish.rnk.dat) <- c("journal","h5ln", "IF13", "IM", "SNIP", "SJR")
mish.avg.rnk <- apply(mish.rnk.dat[,2:6], 1, "mean")
mish.sd.rnk <- apply(mish.rnk.dat[,2:6], 1, "sd")
mish.rnk.dat <- data.frame(mish.rnk.dat, mish.avg.rnk, mish.sd.rnk)
mish.rnk.sort <- mish.rnk.dat[order(mish.rnk.dat[,7],decreasing=F),]
mish.rnk.sort
plot(mish.rnk.sort$mish.avg.rnk, mish.rnk.sort$mish.sd.rnk)

ogyn.h5ln.rnk <- logyn - rank(h5ln.ogyn, ties.method="average") + 1
ogyn.IF13.rnk <- logyn - rank(dat.ogyn$IF2013, ties.method="average") + 1
ogyn.IM.rnk <- logyn - rank(dat.ogyn$IM, ties.method="average") + 1
ogyn.SNIP13.rnk <- logyn - rank(dat.ogyn$SNIP13, ties.method="average") + 1
ogyn.SJR13.rnk <- logyn - rank(dat.ogyn$SJR13, ties.method="average") + 1
ogyn.rnk.dat <- data.frame(dat.ogyn$Journal, ogyn.h5ln.rnk, ogyn.IF13.rnk, ogyn.IM.rnk, ogyn.SNIP13.rnk, ogyn.SJR13.rnk)
colnames(ogyn.rnk.dat) <- c("journal","h5ln", "IF13", "IM", "SNIP", "SJR")
ogyn.avg.rnk <- apply(ogyn.rnk.dat[,2:6], 1, "mean")
ogyn.sd.rnk <- apply(ogyn.rnk.dat[,2:6], 1, "sd")
ogyn.rnk.dat <- data.frame(ogyn.rnk.dat, ogyn.avg.rnk, ogyn.sd.rnk)
ogyn.rnk.sort <- ogyn.rnk.dat[order(ogyn.rnk.dat[,7],decreasing=F),]
ogyn.rnk.sort
plot(ogyn.rnk.sort$ogyn.avg.rnk, ogyn.rnk.sort$ogyn.sd.rnk)


###################
## Bootstrap rank
###################

#dat <- dat.ecol
#dat <- dat.ecol.samp
#dat <- dat.ecol.samp14
#dat <- dat.ecol.samp15
#dat <- subset(dat.ecol.samp15, ecol == 1)
#dat <- subset(dat.ecol.samp15, conserv == 1)
#dat <- subset(dat.ecol.samp14, review == 1)
#dat <- subset(dat.ecol.samp14, biogeo == 1)
#dat <- subset(dat.ecol.samp14, bes == 1)
dat <- subset(dat.ecol.samp15, mar == 1)
#dat <- dat.ener

#dat <- dat.medi
#dat <- dat.mult
#dat <- dat.mish
#dat <- dat.ogyn

ldat <- dim(dat)[1]
no.vec <- 1:ldat
iter <- 10000
itdiv <- iter/100
rnk.mat <- rnk.mat.jk <- matrix(data=NA, nrow=iter, ncol=ldat)
colnames(rnk.mat) <- dat$Journal; colnames(rnk.mat.jk) <- dat$Journal

h5ln.obs <- ldat - rank(dat$h5/log10(dat$n), ties.method="average") + 1
IF.obs <- ldat - rank(dat$IF, ties.method="average") + 1
IM.obs <- ldat - rank(dat$IM, ties.method="average") + 1
SNIP.obs <- ldat - rank(dat$SNIP, ties.method="average") + 1
SJR.obs <- ldat - rank(dat$SJR, ties.method="average") + 1
dat.raw <- data.frame(dat$Journal, h5ln.obs, IF.obs, IM.obs, SNIP.obs, SJR.obs)
avg.rnk <- apply(dat.raw[,2:6], 1, "mean")
dat.obs <- data.frame(dat$Journal, avg.rnk)

for (i in 1:iter) {
  sub.boot <- sample(no.vec, ldat, replace=T)
  dat.boot <- dat[sub.boot,]
  h5ln.rnk <- ldat - rank(dat.boot$h5/log10(dat.boot$n), ties.method="average") + 1
  IF.rnk <- ldat - rank(dat.boot$IF, ties.method="average") + 1
  IM.rnk <- ldat - rank(dat.boot$IM, ties.method="average") + 1
  SNIP.rnk <- ldat - rank(dat.boot$SNIP, ties.method="average") + 1
  SJR.rnk <- ldat - rank(dat.boot$SJR, ties.method="average") + 1
  rnk.boot <- data.frame(dat.boot$Journal, h5ln.rnk, IF.rnk, IM.rnk, SNIP.rnk, SJR.rnk)
  colnames(rnk.boot) <- c("journal","h5ln", "IF", "IM", "SNIP", "SJR")
  boot.avg.rnk <- apply(rnk.boot[,2:6], 1, "mean")
  rnk.boot.dat <- data.frame(rnk.boot, boot.avg.rnk)
  rnk.boot.sort <- rnk.boot.dat[order(rnk.boot.dat[,7],decreasing=F),]
  rnk.boot.sort.unique <- unique(rnk.boot.sort)
  rnk.boot.sort.unique
  
  # store bootstrapped ranks
  lboot <- dim(rnk.boot.sort.unique)[1]
  sub.stor <- rep(0,lboot)
  for (b in 1:lboot) {
    sub.stor[b] <- which(dat$Journal == rnk.boot.sort.unique$journal[b])
  }
  
  rnk.mat[i,sub.stor] <- rnk.boot.sort.unique$boot.avg.rnk
  if (i %% itdiv==0) print(i)
}
rnk.med <- apply(rnk.mat, 2, "median", na.rm=T)
rnk.lo <- apply(rnk.mat, 2, "quantile", na.rm=T, probs=0.025)
rnk.up <- apply(rnk.mat, 2, "quantile", na.rm=T, probs=0.975)

# kappa method
kappa <- 2
rnk.update <- rnk.mat
kappa.n <- 5

for (k in 1:kappa.n) {
  boot.avg <- apply(rnk.update, 2, "mean", na.rm=T)
  boot.sd <- apply(rnk.update, 2, "sd", na.rm=T)
  subset.mat <- rnk.new <- matrix(data=NA, nrow=iter, ncol=ldat)
  
  for (l in 1:ldat) {
    subset.sub <- which(rnk.update[,l] > (boot.avg[l] - (kappa*boot.sd[l])) & rnk.mat[,l] < (boot.avg[l] + (kappa*boot.sd[l])))
    lsubset <- length(subset.sub)
    subset.mat[1:lsubset,l] <- subset.sub
  } 
  for (l in 1:ldat) {
    rnk.new.vals <- rnk.mat[na.omit(subset.mat[,l]),l]
    lnew <- length(rnk.new.vals)
    rnk.new[1:lnew,l] <- rnk.new.vals
  }
  rnk.update <- rnk.new
  print(k)
}

rnk.upd.med <- apply(rnk.update, 2, "median", na.rm=T)
rnk.upd.lo <- apply(rnk.update, 2, "quantile", na.rm=T, probs=0.025)
rnk.upd.up <- apply(rnk.update, 2, "quantile", na.rm=T, probs=0.975)

# jackknife
jnl.se.rnk <- jnl.up.rnk <- jnl.lo.rnk <- rep(0,ldat)
theta <- function(x){mean(x)}

for (i in 1:ldat) {
  jnl.jk <- jackknife(as.numeric(dat.raw[i,2:6]), theta)
  jnl.se.rnk[i] <- jnl.jk$jack.se
  jnl.up.rnk[i] <- quantile(jnl.jk$jack.values, probs=0.975)
  jnl.lo.rnk[i] <- quantile(jnl.jk$jack.values, probs=0.025)
}

out.jk <- data.frame(dat.obs,jnl.lo.rnk,jnl.up.rnk,jnl.se.rnk)
out.jk
out.jk.sort <- out.jk[order(out.jk[,2], decreasing=F),]
#setwd("~/Documents/Papers/Other/Journal Ranks/ms/PLoS One/R1/results/")
#write.table(out.jk.sort,file="ecol.jk.csv",sep=",",dec = ".", row.names = T,col.names = TRUE)

dat.out <- data.frame(dat.obs$avg.rnk,rnk.lo,rnk.med,rnk.up,rnk.upd.lo,rnk.upd.med,rnk.upd.up,jnl.up.rnk,jnl.lo.rnk,jnl.se.rnk)
colnames(dat.out) <- c("avg", "lo", "med", "up","lo.k", "med.k", "up.k", "jk.up", "jk.lo", "jk.se")
dat.sort <- dat.out[order(dat.out[,6],decreasing=F),]
dat.sort
setwd("~/Documents/Papers/Other/Journal Ranks/ms/PLoS One/R1/results/")
write.table(dat.sort,file="mar15.sort.csv",sep=",",dec = ".", row.names = T,col.names = TRUE)

## plot error bars
errbar(rownames(dat.sort)[1:10], dat.sort$med.k[1:10], yplus=dat.sort$up.k[1:10], yminus=dat.sort$lo.k[1:10],xlab="rank", cex.axis = 0.5)
errbar(rownames(dat.sort), dat.sort$med.k, yplus=dat.sort$up.k, yminus=dat.sort$lo.k,xlab="rank", cex.axis = 0.5)

IFrnk.dat <- data.frame(dat$Journal,dat$IF)
colnames(IFrnk.dat) <- c("Journal","IF")
datsort <- data.frame(dat.sort, rownames(dat.sort))
colnames(datsort) <- c("avg", "lo", "med", "up", "lo.k", "med.k", "up.k", "jk.up", "jk.lo", "jk.se", "Journal")
IFcom.merge <- merge(IFrnk.dat,datsort,"Journal")
cor(IFcom.merge[,2], IFcom.merge[,8], method="pearson")
plot(IFcom.merge[,2], IFcom.merge[,8], xlab="IF", ylab="composite rank")



## Compare mean- vs. median-calculated ranks
setwd("~/Documents/Papers/Other/Journal Ranks/ms/PLoS One/R1/results/")
ecol.med <- read.table("ecol.sort.csv", header=T, sep = ",") 
colnames(ecol.med) <- c("journal","avg.md","lo.md","med.md","up.md","lo.k.md","med.k.md","up.k.md","jk.up.md","jk.lo.md","jk.se.md")
setwd("~/Documents/Papers/Other/Journal Ranks/results/")
ecol.mn <- read.table("ecol.sort.csv", header=T, sep = ",") 
ecol.comp <- merge(ecol.mn, ecol.med, "journal")
setwd("~/Documents/Papers/Other/Journal Ranks/ms/PLoS One/R1/results/")
write.table(ecol.comp,file="ecol.comp.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)

setwd("~/Documents/Papers/Other/Journal Ranks/ms/PLoS One/R1/results/")
medi.med <- read.table("medi.sort.csv", header=T, sep = ",") 
colnames(medi.med) <- c("journal","avg.md","lo.md","med.md","up.md","lo.k.md","med.k.md","up.k.md","jk.up.md","jk.lo.md","jk.se.md")
setwd("~/Documents/Papers/Other/Journal Ranks/results/")
medi.mn <- read.table("medi.samp.sort.csv", header=T, sep = ",") 
medi.comp <- merge(medi.mn, medi.med, "journal")
setwd("~/Documents/Papers/Other/Journal Ranks/ms/PLoS One/R1/results/")
write.table(medi.comp,file="medi.comp.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)

setwd("~/Documents/Papers/Other/Journal Ranks/ms/PLoS One/R1/results/")
mult.med <- read.table("mult.sort.csv", header=T, sep = ",") 
colnames(mult.med) <- c("journal","avg.md","lo.md","med.md","up.md","lo.k.md","med.k.md","up.k.md","jk.up.md","jk.lo.md","jk.se.md")
setwd("~/Documents/Papers/Other/Journal Ranks/results/")
mult.mn <- read.table("mult.samp.sort.csv", header=T, sep = ",") 
mult.comp <- merge(mult.mn, mult.med, "journal")
setwd("~/Documents/Papers/Other/Journal Ranks/ms/PLoS One/R1/results/")
write.table(mult.comp,file="mult.comp.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)



## Cluster Analyses
#library(OpenRepGrid)
library(pvclust)
dat <- dat.ecol.samp[,c(3,5,7,9:10,12)]
#dat <- dat.ecol[1:30,c(3,5,7,9:10,12)]
#dat <- dat.medi[1:30,c(3,5,7,9:10,12)]
#dat <- dat.mult[1:30,c(3,5,7,9:10,12)]

## h5/n
h5n <- dat[,2]/log10(dat[,1])
dat <- data.frame(dat,h5n)
dat <- dat[,-c(1,2)]
dat

rownames(dat) <- dat.ecol.samp$Journal
#rownames(dat) <- dat.ecol$Journal[1:30]
#rownames(dat) <- dat.medi$Journal[1:30]
#rownames(dat) <- dat.mult$Journal[1:30]

out.agnes <- agnes(dat, diss=F, metric="euclidian", stand=T, method="complete", trace.lev=2)
plot(out.agnes, main="")

out.diana <- diana(dat, diss=F, metric="euclidian", stand=T, trace.lev=2)
plot(out.diana, main="")

# scale data
dat.scale <- scale(dat, center=T, scale=T)
dat.tscale <- as.data.frame(t(dat.scale))
dat.dist <- dist(dat.scale)
#dat.daisy <- daisy(dat, metric="euclidian", stand=T)
out.hclust <- hclust(dat.dist, method="complete")
plot(out.hclust,labels=rownames(dat),main="",sub="")

out.pvclust <- pvclust(dat.tscale, method="complete", nboot=10000, method.dist="euclidean")
plot(out.pvclust)
pvrect(out.pvclust)
out.pp <- pvpick(out.pvclust)
out.pp

## PCA
#library(rda)
library(Rcmdr)
library(vegan)
z <- rda(dat.scale, data = dat.scale, scale=T, na.action="na.omit")
plot(z, type="text")

## summaries
summary(z)
summ.out <- summary(z)
pca.out <- summ.out$sites
#setwd("~/Documents/Papers/Other/Journal Ranks/results/")
#write.table(pca.out,file="pca.out.csv",sep=",",dec = ".", row.names = T,col.names = TRUE)





## Survey
dat.survey <- read.table("survey.csv", header=T, sep=",")
dat.surv <- subset(dat.survey, keep==1)
dat.surv <- subset(dat.surv, articles == "51-100" | articles == "> 100") # choose only people with > 50 pubs

col.vec <- colnames(dat.surv)
journal.vec <- colnames(dat.surv)[3:27]
ljourn <- length(journal.vec)

for (j in 1:ljourn) {
  sub.j <- which(colnames(dat.surv) == journal.vec[j])
  choice1 <- gsub("1-Elite", 1, dat.surv[,sub.j])
  choice2 <- gsub("2-Prestigious", 2, choice1)
  choice3 <- gsub("3-Reputable", 3, choice2)
  choice4 <- gsub("4-Respectable", 4, choice3)  
  dat.surv[,sub.j] <- as.numeric(gsub("5-Other", 5, choice4))
}

choice.mean <- apply(dat.surv[,3:27],2,"mean")
choice.sd <- apply(dat.surv[,3:27],2,"sd")
choice.lo <- choice.mean - choice.sd
choice.up <- choice.mean + choice.sd
jrank <- rank(choice.mean, ties.method="average")
choice.out <- data.frame(jrank,choice.mean,choice.sd,choice.lo,choice.up)
colnames(choice.out) <- c("rank", "mean", "sd", "lo", "up")
choice.sort <- choice.out[order(choice.out[,1], decreasing=F),]
choice.sort
errbar(rownames(choice.sort)[1:15], choice.sort$mean[1:15], yplus=choice.sort$up[1:15], yminus=choice.sort$lo[1:15],xlab="rank", cex.axis = 0.5)
#setwd("~/Documents/Papers/Other/Journal Ranks/results/")
#write.table(choice.sort,file="choice.sort.m50.csv",sep=",",dec = ".", row.names = T,col.names = TRUE)

## Correlation between ecol.samp & choice.sort
e.choice <- read.table("choice.sort.m50.csv", header=T, sep=",")
e.samp <- read.table("ecol.samp.kappa.sort.csv", header=T, sep=",")
e.samp$journal <- e.samp$X

# merge
e.merge <- merge(e.samp, e.choice, by="journal")

# iterate (sample within range for Spearman's correlation)
iter <- 1000

xs <- data.frame(e.merge$lo.k.x,e.merge$up.k.x)
colnames(xs) <- c("lo","up")
ys <- data.frame(e.merge$lo.k.y,e.merge$up.k.y)
colnames(ys) <- c("lo","up")
xs.mat <- matrix(data=0, nrow=iter, ncol=dim(xs)[1])
ys.mat <- matrix(data=0, nrow=iter, ncol=dim(ys)[1])

for (r in 1:dim(xs)[1]) {
    xs.mat[,r] <- runif(iter, xs[r,1], xs[r,2])
    ys.mat[,r] <- runif(iter, ys[r,1], ys[r,2])
}

rho.vec <- rep(0,iter)
for (i in 1:iter) {
  rho.vec[i] <- cor(xs.mat[i,], ys.mat[i,], method="spearman")
}
rho.med <- median(rho.vec)
rho.lo <- quantile(rho.vec,probs=0.025)
rho.up <- quantile(rho.vec,probs=0.975)
print(c(round(rho.lo,3),round(rho.med,3),round(rho.up,3)))



## Bootstrap
dat <- dat.surv
ldat <- dim(dat)[1]
no.vec <- 1:ldat
iter <- 10000
itdiv <- iter/100
rnk.mat <- matrix(data=NA, nrow=iter, ncol=ljourn)
colnames(rnk.mat) <- journals.vec

for (i in 1:iter) {
  sub.boot <- sample(no.vec, ldat, replace=T)
  dat.boot <- dat[sub.boot,]
  mean <- apply(dat.boot[,3:27],2,"mean")
  rank <- rank(mean, ties.method="average")
  rnk.boot.dat <- data.frame(mean, rank)
  rnk.mat[i,] <- rnk.boot.dat$rank
  if (i %% itdiv==0) print(i)
}
rnk.med <- apply(rnk.mat, 2, "median", na.rm=T)
rnk.lo <- apply(rnk.mat, 2, "quantile", na.rm=T, probs=0.025)
rnk.up <- apply(rnk.mat, 2, "quantile", na.rm=T, probs=0.975)


# kappa method
kappa <- 2
kappa.n <- 5
rnk.update <- rnk.mat

for (k in 1:kappa.n) {
  boot.avg <- apply(rnk.update, 2, "mean", na.rm=T)
  boot.sd <- apply(rnk.update, 2, "sd", na.rm=T)
  sd.zero <- which(boot.sd == 0)
  if (length(sd.zero) > 0) {
    boot.sd[sd.zero] <- as.numeric(min(boot.sd[-sd.zero]))
  }
  subset.mat <- rnk.new <- matrix(data=NA, nrow=iter, ncol=ljourn)
  
  for (l in 1:ljourn) {
    subset.sub <- which(rnk.update[,l] > (boot.avg[l] - (kappa*boot.sd[l])) & rnk.mat[,l] < (boot.avg[l] + (kappa*boot.sd[l])))
    lsubset <- length(subset.sub)
    subset.mat[1:lsubset,l] <- subset.sub
  } 
  for (l in 1:ljourn) {
    rnk.new.vals <- rnk.mat[na.omit(subset.mat[,l]),l]
    lnew <- length(rnk.new.vals)
    rnk.new[1:lnew,l] <- rnk.new.vals
  }
  rnk.update <- rnk.new
  print(k)
}

rnk.upd.med <- apply(rnk.update, 2, "median", na.rm=T)
rnk.upd.lo <- apply(rnk.update, 2, "quantile", na.rm=T, probs=0.025)
rnk.upd.up <- apply(rnk.update, 2, "quantile", na.rm=T, probs=0.975)

dat.out <- data.frame(choice.mean,jrank,rnk.upd.lo,rnk.upd.med,rnk.upd.up)
colnames(dat.out) <- c("mean", "mn.rnk", "lo.k", "med.k", "up.k")
dat.sort <- dat.out[order(dat.out[,4],decreasing=F),]
dat.sort
#setwd("~/Documents/Papers/Other/Journal Ranks/results/")
#write.table(dat.sort,file="surv.sort.kappa.m50.csv",sep=",",dec = ".", row.names = T,col.names = TRUE)

## plot error bars
errbar(rownames(dat.sort)[1:20], dat.sort$med.k[1:20], yplus=dat.sort$up.k[1:20], yminus=dat.sort$lo.k[1:20],xlab="rank", cex.axis = 0.5)

