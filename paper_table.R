options(scipen=9)
source('function_orbit_v6.R')
library(xtable)
load('results/clone5pc.Robj')
tstat <- sapply(1:nrow(Tenc), function(i) clone.stat(Tenc[i,]))
dstat <- sapply(1:nrow(Denc), function(i) clone.stat(Denc[i,]))
vstat <- sapply(1:nrow(Venc), function(i) clone.stat(Venc[i,]))
xenc.stat <- sapply(1:nrow(Xenc), function(i) clone.stat(Xenc[i,]))
yenc.stat <- sapply(1:nrow(Yenc), function(i) clone.stat(Yenc[i,]))
zenc.stat <- sapply(1:nrow(Zenc), function(i) clone.stat(Zenc[i,]))
tab <- read.csv('results/enc5pc.csv')
tout <- cbind(Tenc[,1],t(tstat[c('x5per','x95per'),]))
colnames(tout)[1] <- 'tenc'
d5per <- sapply(1:nrow(Denc), function(i) min(Denc[i,]))
#d95per <- sqrt(xenc.stat['x95per',]^2+yenc.stat['x95per',]^2+zenc.stat['x95per',]^2)
d95per <- dstat['x95per',]
#dout <- cbind(Denc[,1],t(dstat[c('x5per','x95per'),]))
dout <- cbind(Denc[,1],d5per,d95per)
colnames(dout)[1] <- 'denc'
vout <- cbind(Venc[,1],t(vstat[c('x5per','x95per'),]))
colnames(vout)[1] <- 'venc'
inds <- which(tab[,'denc']<2 | tab[,'venc']<10)
val <- cbind(tab[,c('TGAS.ID','HIP.ID','cat')],tout,dout,vout)
val.all <- cbind(tab[,-(1:5)],tout,dout,vout)

####save data
colnames(val) <- c('TGAS','HIP','cat','tenc','tenc.5per','tenc.95per','denc','denc.5per','denc.95per','venc','venc.5per','venc.95per')
write.csv(val,file='enc5pc_final.csv',quote=FALSE,row.names=FALSE)
colnames(val.all) <- c(colnames(tab)[-(1:5)],'tenc','tenc.5per','tenc.95per','denc','denc.5per','denc.95per','venc','venc.5per','venc.95per')
write.csv(val.all,file='enc5pc_all.csv',quote=FALSE,row.names=FALSE)

#####plot paper table
out <- val[inds,]
kk <- sort(out[,'denc'],index.return=TRUE)$ix
names <- c('HIP 21553','HIP 71681','HIP 70890','HIP 71683','HIP 17288','HIP 104539','TGAS 7582-1449-1','HIP 101180','HIP 24608','HIP 86916','HIP 87937','TYC 5855-2215-1','HIP 103749','HIP 113020','HIP 107556','HIP 37766','HIP 51966')
tmp <- cbind(names,out[kk,3:ncol(out)])
pp <- xtable(tmp)
#print(pp,include.rownames=FALSE)
#####astrometric info
#print(xtable(astrometry,digits=c(rep(0,nd),rep(2,12))),include.rownames=FALSE)
source('plot_denc_venc.R')

