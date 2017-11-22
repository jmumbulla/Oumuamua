args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    Popt0 <- as.numeric(args[1])
    type <- args[2]
}else{
    Popt0 <- 8.26
    type <- 'red'
}
#fpdf <- paste0('phase_data_',type,format(Popt*24,digit=3),'h.pdf')
fpdf <- paste0('phase_data_all.pdf')
cat('output pdf:\n',fpdf,'\n')
pdf(fpdf,10,8)
size <- 1.2
#par(mfrow=c(2,2),mar=c(4.5,4.5,1,1),cex.lab=size,cex.axis=size,cex=size)
par(mfrow=c(2,2),mar=c(1,1,4.8,1),cex.lab=size,cex.axis=size,cex=size)
Popts <- c(1,8.10,8.14,8.26)
for(kk in 1:length(Popts)){
xaxt <- 's'
yaxt <- 's'
if(kk==1){
par(mar=c(0,4.8,4.8,0))
xaxt <- 'n'
}else if(kk==3){
par(mar=c(4.8,4.8,0,0))
}else if(kk==4){
par(mar=c(4.8,0,0,4.8))
yaxt <- 'n'
}else if(kk==2){
par(mar=c(0,0,4.8,4.8))
xaxt <- 'n'
yaxt <- 'n'
}
Popt0 <- Popts[kk]
if(Popt0==1){
    load('../../dwarfs/output/U1/keppure_priore0_poly10_Ndata105_quantifyTRUE_1per1_Nw1_U1_all_abs_ind0_1planet_ARMA01_Nsamp1000000_tem1_acc3.2_pretem0.512P0.1d_negLmax17.6.Robj')
    load('../../dwarfs/output/U1/keppure_priore0_poly10_Ndata105_quantifyTRUE_1per1_Nw1_U1_all_abs_ind0_1planet_ARMA01_Nsamp1000000_tem1_acc3.2_pretem0.512P0.1d_negLmax17.6_optpar.Robj')
#    Popt <- par.stat['mean',1]
}else if(Popt0==2){
    load('../../dwarfs/output/U1/keppure_priore0_poly10_Ndata105_quantifyTRUE_1per1_Nw1_U1_all_abs_ind0_1planet_ARMA01_Nsamp1000000_tem1_acc9.1_pretem1P0.4d_negLmax20.4.Robj')
}else{
    Popt <- Popt0/24#d
}
if(Popt>6/24){
    ofac <- 1
}else{
    ofac <- 2
}
tab <- read.table('U1_all.txt',header=TRUE)
out <- read.table('U1_all_abs.txt',header=TRUE)
#Popt <- 8.26/24
if(Popt0==1 | Popt0==2){
    a <- par.opt.post[grep(paste0('a\\d_',k1),names(par.opt.post))]
    b <- par.opt.post[grep(paste0('\\db',k1),names(par.opt.post))]
    trend <- cal.trend(a=a,b=b,t=(trv-min(trv))/time.unit)
    if(type=='red'){
        mags <- RV2+b
    }else{
        mags <- out[,2]-a*(trv-min(trv))/time.unit
#        mags <- out[,2]
    }
}else{
    x <- out[,1]
    y <- out[,2]
#    mags <- out[,2]
#    val <- lm(y ~ poly(x,2))
#    mags <- residuals(val)+as.numeric(val$coefficients[1])
}
#Popt <- 8.26/24
#x <- out[,1]
#y <- out[,2]
#emag <- out[,3]
emag <- out[,3]
#mags <- out[,2]
#val <- lm(y ~ poly(x,1))
#mags <- residuals(val)#residual
#par(mfrow=c(2,2))
##polynomial fit
#plot(out[,1],out[,2],xlab='MJD',ylab='Rmag')
#tsim <- seq(min(out[,1]),max(out[,1]),length.out=100)
#lines(tsim,predict(val,data.frame(x=tsim)),col='red')
##periodic function fit
t <- out[,1]
t <- (t%%(ofac*Popt))*24#h
ns <- c('APO','NOT','WIYN','DCT','GNT','WHT')
cs <- c('black','brown','red','blue','green','orange')
for(k in 1:length(ns)){
    ind <- which(tab[,'telescope']==ns[k])
    if(k==1){
        plot(t[ind],mags[ind],xlab='',ylab='',pch=20,ylim=rev(range(mags,24,20.5)),xlim=range(t),col=cs[k],xaxt=xaxt,yaxt=yaxt)
    }else{
        points(t[ind],mags[ind],col=cs[k],pch=20)
    }
    arrows(t[ind],mags[ind]-emag[ind],t[ind],mags[ind]+emag[ind],length=0.03,angle=90,code=3,col=cs[k])
}
if(Popt0==1|Popt0==2){
    ref <- 'This work'
}else if(Popt0==8.26){
    ref <- 'Jewitt et al. 2017'
}else if(Popt0==8.10){
    ref <- 'Bannister et al. 2017'
}else if(Popt0==8.14){
    ref <- 'Bolin et al. 2014'
}
if(kk==2){
legend('topright',xpd=NA,inset=c(-0.25,0),legend=ns,col=cs,pch=20)
}
if(kk==1){
mtext(side=2,outer=TRUE,line=-2,text='Absolute Red Magnitude',cex=1.5)
}
if(kk==3){
mtext(side=1,outer=TRUE,line=-2,text='Time [hours]',cex=1.5)
}
legend('topright',legend=paste0(format(Popt*ofac*24,digit=3),' hour (',ref,')'),bty='n')
}
dev.off()
