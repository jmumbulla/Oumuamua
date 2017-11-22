tab1 <- read.table('NOT_WIYN.txt',header=TRUE)
tab2 <- read.table('APO_DCT.txt',header=TRUE)
tab3 <- read.table('OSSOS.txt',header=TRUE)
#name <- 'U1'
name <- 'NW'
#name <- 'AD'
#name <- 'OS'
ns <- c('NW','AD','OS','NW_AD','NW_OS','AD_OS','U1')
###ref. Jewitt et al. 2017
bv <- 0.7
ebv <- 0.06
vr <- 0.45
evr <- 0.05
br <- 1.15
ebr <- 0.05
##convert to R mag
tab1[,'date'] <- tab1[,'date']-25+58051.5
colnames(tab1)[1] <- 'MJD'
magR <- tab1[,'mag']
indB <- which(tab1[,'filter']=='B')
indV <- which(tab1[,'filter']=='V')
magR[indB] <- tab1[indB,'mag']-br
magR[indV] <- tab1[indV,'mag']-vr
##add instruments for APO and DCT data
ins <- rep('DCT',nrow(tab2))
ind <- which(tab2[,1]<58055.6)
ins[ind] <- 'APO'
##process the OSSOS data
gr.GNT <- 0.47
gr.WHT <- 0.63
ri.WHT <- 0.36
rj.GNT <- 1.2
magr <- tab3[,'mag']
ind.gg <- which(tab3[,'filter']=='gg')
magr[ind.gg] <- tab3[ind.gg,'mag']-gr.GNT
ind.j <- which(tab3[,'filter']=='j')
magr[ind.j] <- tab3[ind.j,'mag']+rj.GNT
ind.gi <- which(tab3[,'filter']=='gi')
magr[ind.gi] <- tab3[ind.gi,'mag']-gr.WHT
ind.ii <- which(tab3[,'filter']=='ii')
magr[ind.ii] <- tab3[ind.ii,'mag']+ri.WHT
##combine too catalogs
out1 <- cbind(tab1[,'MJD'],magR,tab1[,c('emag','telescope')])
out2 <- cbind(tab2[,'MJD'],tab2[,c('magR','emagR')],ins)
out3 <- cbind(tab3[,'epoch'],magr,tab3[,c('emag','telescope')])
colnames(out3) <- colnames(out2) <- colnames(out1)
#out <- rbind(out1,out2,out3)
for(name in ns){
if(name=='NW'){
out <- out1
}else if(name=='AD'){
out <- out2
}else if(name=='OS'){
out <- out3
}else if(name=='NW_AD'){
out <- rbind(out1,out2)
}else if(name=='NW_OS'){
out <- rbind(out1,out3)
}else if(name=='AD_OS'){
out <- rbind(out2,out3)
}else if(name=='U1'){
out <- rbind(out1,out2,out3)
}
colnames(out) <- c('MJD','magR','emagR','telescope')
index <- sort(as.numeric(out[,'MJD']),index.return=TRUE)$ix
out <- out[index,]
out[,1] <- as.numeric(out[,1])
f1 <- paste0(name,'_all.txt')
f2 <- paste0(name,'_all_simple.txt')
write.table(out,file=f1,row.names=FALSE,quote=FALSE)
write.table(out[,1:3],file=f2,row.names=FALSE,quote=FALSE)
###subtract the trend
x <- out[,1]
y <- out[,2]
val <- lm(y ~ poly(x,3))
mags <- residuals(val)#residual
fpdf <- paste0(name,'_data.pdf')
cat('output pdf:\n',fpdf,'\n')
pdf(fpdf,8,8)
par(mfrow=c(2,2))
##polynomial fit
plot(out[,1],out[,2],xlab='MJD',ylab='Rmag')
tsim <- seq(min(out[,1]),max(out[,1]),length.out=100)
lines(tsim,predict(val,data.frame(x=tsim)),col='red')
##periodic function fit
for(j in 1:2){
mag <- mags
if(j==1){
    t <- out[,1]
    plot(t,mag,xlab='MJD',ylab='Rmag')
}else{
    t <- (t%%0.1455)*24
    plot(t,mag,xlab='phase-folded time [h]',ylab='Rmag')
}
ind.not <- which(out[,'telescope']=='NOT')
ind.wiyn <- which(out[,'telescope']=='WIYN')
ind.dct <- which(out[,'telescope']=='DCT')
ind.gnt <- which(out[,'telescope']=='GNT')
ind.wht <- which(out[,'telescope']=='WHT')
points(t[ind.not],mag[ind.not],col='yellow')
points(t[ind.wiyn],mag[ind.wiyn],col='red')
points(t[ind.dct],mag[ind.dct],col='blue')
points(t[ind.gnt],mag[ind.gnt],col='green')
points(t[ind.wht],mag[ind.wht],col='cyan')
}
dev.off()
tmp <- cbind(out[,1],mags,out[,3])
colnames(tmp) <- c('MJD','magR','emagR')
f3 <- paste0(name,'_all_residual.txt')
cat('output data:\n',f3,'\n\n')
write.table(tmp,file=f3,row.names=FALSE,quote=FALSE)
}
