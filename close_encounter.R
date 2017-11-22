refine <- TRUE
ds <- c()
ts <- c()
dmin <- c()
tmin <- c()
for(j in 1:Nt){
    xx <- xs[,j,1]-xs[,1,1]
    yy <- ys[,j,1]-ys[,1,1]
    zz <- zs[,j,1]-zs[,1,1]
    vxx <- vxs[,j,1]-vxs[,1,1]
    vyy <- vys[,j,1]-vys[,1,1]
    vzz <- vzs[,j,1]-vzs[,1,1]
    dd <- sqrt(xx^2+yy^2+zz^2)*1e3
    ind <- which.min(dd)
    if(refine){
        tmp <- lp(times[ind],r=c(xx[ind],yy[ind],zz[ind]),v=c(vxx[ind],vyy[ind],vzz[ind]))
        if(j==1){
            ts <- cbind(ts,c(0,times))
            ds <- cbind(ds,c(0,dd))
        }else{
            tt <- c(times,tmp[2])
            indt <- sort(tt,index.return=TRUE)$ix
            ts <- cbind(ts,tt[indt])
            ds <- cbind(ds,c(dd,tmp[1]*1e3)[indt])
        }
        dmin <- c(dmin,tmp[1]*1e3)
        tmin <- c(tmin,tmp[2])
    }else{
        ts <- cbind(ts,times)
        ds <- cbind(ds,dd)
        dmin <- c(dmin,dd[ind])
        tmin <- c(tmin,times[ind])
    }
    if(min(dd)<1){
        cat('denc=',dd[ind],'pc\n')
        cat('tenc=',times[ind],'Myr\n')
        cat('ind.star=',ind.star[j],'\n\n')
    }
}
#names <- c('Sun','HIP 74273','HIP 54035','HIP 72511','HIP 20074','J064422.74+415517.7','HIP 104965','HIP 24284','HIP 35296')
names <- c('Sun','38 G. Circini','HD 95735','LHS 380','HD 27553','J064422.74+415517.7','HD 202027','HIP 24284','HD 57095')
pdf('encounter_trajectory.pdf',8,4)
par(mar=c(4.2,4.2,1,17))
dproxima <- 13/206.265#pc
cc <- rainbow(Nt-2, alpha = 1)
plot(ts[,3],ds[,3],xlab='Time [Myr]',ylab='Periastron [pc]',col=cc[1],ylim=c(0.01,1.5),xlim=c(-3,1),type='l')#log='y'
abline(h=1,lty=2)
#text(tmin[3],dmin[3],labels='Sun',col=cc[1])
#abline(h=,lty=3)
x3 <- seq(-4,2,by=0.1)
y3 <- seq(-1,dproxima,by=0.1)
Ny <- length(y3)
Nx <- length(x3)
polygon(x=c(x3,rep(2,Ny),rev(x3),rep(-4,Ny)), y =c(rep(-1,Nx),y3,rep(dproxima,Nx),rev(y3)) ,col=tcol('grey',60),border=NA)
text(-1,dproxima/2,labels='Orbit of Proxima',cex=0.8)
for(j in 4:Nt){
    lines(ts[,j],ds[,j],col=cc[j-2])
}
signs <- rep('+',Nt)
signs[tmin<0] <- '-'
legend('topright',inset=c(-0.9,0),legend=paste0(names,' (',signs[3:Nt],format(abs(tmin[3:Nt]),digit=1),' Myr)'),col=cc,xpd=NA,bty='n',lty=rep(1,length(cc)))
dev.off()

