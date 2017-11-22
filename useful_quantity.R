###In heliocentric frame 
if(TRUE){
    cat('K=',K[1],'\n')
    cat('V=',V[1],'\n')
    cat('vesc=',vesc[1]*auyr2kms,'km/s\n')
    Ndig <- 3
    Xgal <- xyzuvw[1,(1:Nt)*6-5]
    Ygal <- xyzuvw[1,(1:Nt)*6-4]
    Zgal <- xyzuvw[1,(1:Nt)*6-3]
    cat('Xgal=',format(Xgal,digit=Ndig),'pc\n')
    cat('Ygal=',format(Ygal,digit=Ndig),'pc\n')
    cat('Zgal=',format(Zgal,digit=Ndig),'pc\n')
    DCom <- r/206265#pc
    cat('Dcom=',format(DCom[1,],digit=Ndig),'pc\n')
    cat('dDcom=',format(apply(DCom,2,'sd'),digit=Ndig),'pc\n')
    U <- xyzuvw[,(1:Nt)*6-2,drop=FALSE]
    V <- xyzuvw[,(1:Nt)*6-1,drop=FALSE]
    W <- xyzuvw[,(1:Nt)*6,drop=FALSE]
    cat('U=',format(U[1,],digit=Ndig),'km/s\n')
    cat('dU=',format(apply(U,2,'sd'),digit=Ndig),'km/s\n') 
    cat('V=',format(V[1,],digit=Ndig),'km/s\n')
    cat('dV=',format(apply(V,2,'sd'),digit=Ndig),'km/s\n')
    cat('W=',format(W[1,],digit=Ndig+1),'km/s\n')
    cat('dW=',format(apply(W,2,'sd'),digit=Ndig),'km/s\n')
    cat('DS=',format(dv[1,]*auyr2kms,digit=Ndig),'km/s\n')
    cat('dDS=',format(apply(dv,2,'sd')*auyr2kms),'km/s\n')
    pm <- hg2pm(xyzuvw[1,(1:Nt)*6-2]/auyr2kms,xyzuvw[1,(1:Nt)*6-1]/auyr2kms,xyzuvw[1,(1:Nt)*6]/auyr2kms,obs.all[,'ra']*pi/180,obs.all[,'dec']*pi/180,1/obs.all[,'plx'])
    cat('pmra [mas/yr]:',pm[,'pmra'],'\n')
    cat('pmdec [mas/yr]:',pm[,'pmde'],'\n')
    cat('obs pmra [mas/yr]:',obs.all[,'pmra'],'\n')
    cat('obs pmdec [mas/yr]:',obs.all[,'pmdec'],'\n')
    if(FALSE){
        cat('\n ',data.type,' results:\n')
        xyzuvw2 <- xyzuvw
        xyzuvw2[,4] <- -5.71
        xyzuvw2[,5] <- -8.26 
        xyzuvw2[,6] <- -11.04
        pm2 <- hg2pm(xyzuvw2[1,(1:Nt)*6-2]/auyr2kms,xyzuvw2[1,(1:Nt)*6-1]/auyr2kms,xyzuvw2[1,(1:Nt)*6]/auyr2kms,obs.all[,'ra']*pi/180,obs.all[,'dec']*pi/180,1/obs.all[,'plx'])
        cat('pmra [mas/yr]:',pm2[,'pmra'],'\n')
        cat('pmdec [mas/yr]:',pm2[,'pmde'],'\n')
    }
}
####plot
fname <- paste0('testAB_',data.type,'_data.pdf')
cat('output pdf:\n',fname,'\n')
pdf(fname,4,4)
size <- 1.0
par(mar=c(4.2,4.2,1,1),cex.lab=size,cex.axis=size,cex=size)
dvb <- sqrt(dvx^2+dvy^2+dvz^2)
hist(Energy,xlab=expression("Orbital energy ["*M[s]*au^2*yr^{-2}*"]"),ylab='Frequency',breaks=500,xlim=c(min(Energy),0.02),main='')
abline(v=0,col='red')
dev.off()
break()
