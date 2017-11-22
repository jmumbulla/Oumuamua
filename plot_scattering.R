####important lesson: For nearby stars, the proper motion, pmra and pmdec, depend on the location of the star. For each time step, the change of position will change the values of pmra and pmdec and distribute the total velocity into three new directions, even though the motions are linear. Therefore, it is not surprising to see that the proper motions of nearby stars in the helio-static frame are different from the ones assuming zero radial velocity and keeping the Sun static. 
NN <- round(tmax-tmin)
#folder <- 'results'
folder <- '/smp3/ffeng/scattering/Oumuamua_result'
if(!file.exists(folder)) folder <- 'results'
pdf.name <- paste0(folder,'/I1_orbit',sim.enc,'_tmin',tmin,'_tmax',tmax,'_dt',dt,'ind',ind11,'_',ind22,'_LA',LA,'.pdf')
cat('output pdf:\n',pdf.name,'\n')
if(!LA){
pdf(pdf.name,12,12)
par(mfrow=c(3,3))
#times <- out[,1,1]
x <- (out[,6*(2:Nt-1)+2,1]-out[,2,1])*1e3#pc
y <- (out[,6*(2:Nt-1)+3,1]-out[,3,1])*1e3#pc
z <- (out[,6*(2:Nt-1)+4,1]-out[,4,1])*1e3#pc
vx <- (out[,6*(2:Nt-1)+5,1]-out[,5,1])*kpcmyr2kms#km/s
vy <- (out[,6*(2:Nt-1)+6,1]-out[,6,1])*kpcmyr2kms#km/s
vz <- (out[,6*(2:Nt-1)+7,1]-out[,7,1])*kpcmyr2kms#km/s
R <- sqrt(x^2+y^2)
v <- sqrt(vx^2+vy^2+vz^2)
if(is.null(dim(R))){
    R <- matrix(R)
    v <- matrix(v)
    x <- matrix(x)
    y <- matrix(y)
    z <- matrix(z)
    vx <- matrix(vx)
    vy <- matrix(vy)
    vz <- matrix(vz)
}
for(j in 1){
    tit <- paste0('Nt=',j)
    plot(times,R[,j],xlab='T [Myr]',ylab='R [pc]',type='l',main=tit)
    plot(times,v[,j],xlab='T [Myr]',ylab='v [km/s]',type='l',main=tit)
    plot(out[,8,1]*1e3,out[,9,1]*1e3,xlab='X [pc]',ylab='Y [pc]',type='l',main=tit)
    plot(out[,8,1]*1e3,out[,10,1]*1e3,xlab='X [pc]',ylab='Z [pc]',type='l',main=tit)
    plot(out[,11,1]*kpcmyr2kms,out[,12,1]*kpcmyr2kms,xlab='Vx [km/s]',ylab='Vy [km/s]',type='l',main=tit)
    plot(out[,11,1]*kpcmyr2kms,out[,13,1]*kpcmyr2kms,xlab='Vx [km/s]',ylab='Vz [km/s]',type='l',main=tit)
    plot(x[,j],y[,j],xlab='x [pc]',ylab='y [pc]',type='l',main=tit)
    plot(x[,j],z[,j],xlab='x [pc]',ylab='z [pc]',type='l',main=tit)
    plot(y[,j],z[,j],xlab='y [pc]',ylab='z [pc]',type='l',main=tit)
    plot(vx[,j],vy[,j],xlab='vx [km/s]',ylab='vy [km/s]',type='l',main=tit)
    plot(vy[,j],vz[,j],xlab='vy [km/s]',ylab='vz [km/s]',type='l',main=tit)
    plot(vx[,j],vz[,j],xlab='vx [km/s]',ylab='vz [km/s]',type='l',main=tit)
    ##bl
    bl <- xyz2bl.vec(x,y,z)
    bl.v <- xyz2bl.vec(vx,vy,vz)
    ad <- gal2equ(bl[,1],bl[,2])
    ad.v <- gal2equ(bl.v[,1],bl.v[,2])
    plot(bl[,2]*180/pi,bl[,1]*180/pi,xlab='L [deg]',ylab='B [deg]',type='l',main=tit)
    points(bl[nrow(bl),2]*180/pi,bl[nrow(bl),1]*180/pi,col='red')
    plot(ad[,2]*180/pi,ad[,1]*180/pi,xlab='RA [deg]',ylab='DEC [deg]',type='l',main=tit)
    points(ad[nrow(ad),2]*180/pi,ad[nrow(ad),1]*180/pi,col='red')
}
if(N1>0){
xe <- (out[,6*(3:Nt-1)+2,1]-out[,6*(2-1)+2,1])*1e3#pc
ye <- (out[,6*(3:Nt-1)+3,1]-out[,6*(2-1)+3,1])*1e3#pc
ze <- (out[,6*(3:Nt-1)+4,1]-out[,6*(2-1)+4,1])*1e3#pc
vxe <- (out[,6*(3:Nt-1)+5,1]-out[,6*(2-1)+5,1])*kpcmyr2kms#km/s
vye <- (out[,6*(3:Nt-1)+6,1]-out[,6*(2-1)+6,1])*kpcmyr2kms#km/s
vze <- (out[,6*(3:Nt-1)+7,1]-out[,6*(2-1)+7,1])*kpcmyr2kms#km/s
d <- sqrt(xe^2+ye^2+ze^2)
v <- sqrt(vxe^2+vye^2+vze^2)
if(is.null(dim(d))){
    d <- matrix(d)
    v <- matrix(v)
}
val <- c()
for(kk in 1:ncol(d)){
    plot(times,d[,kk],xlab='T [Myr]',ylab='d [pc]',main=paste0('encounter',kk))
#    plot(times,v[,kk],xlab='T [Myr]',ylab='v [km/s]',main=paste0('encounter',kk))
    ind.min <- which.min(d[,kk])
    cat('id=',out.all[kk,'id'],'\n')
    cat('denc=',d[ind.min,kk],'\n')
    cat('tenc=',times[ind.min],'\n')
    cat('venc=',v[ind.min],'\n')
    val <- rbind(val,c(d[ind.min,kk],times[ind.min],v[ind.min]))
}
colnames(val) <- c('denc','tenc','venc')
out.final <- cbind(val,out.all[ind11:ind22,])
save(out.final,file=gsub('.pdf','.Robj',pdf.name))
}
#plot(bl.v[,2]*180/pi,bl.v[,1]*180/pi,xlab='Lv [deg]',ylab='Bv [deg]',type='l')
#points(bl.v[nrow(bl.v),2]*180/pi,bl.v[nrow(bl.v),1]*180/pi,col='red')
dev.off()
}else{
    pdf(pdf.name,12,8)
    par(mfrow=c(2,3))
    plot(out[,2,1],out[,3,1],xlab='X [kpc]',ylab='Y [kpc]',type='l')
    uvw <- calc.vp(out[,2,1],out[,3,1],out[,4,1],out[,5,1],out[,6,1],out[,7,1],parameters=par.gal)*kpcmyr2kms
    plot(out[,1,1],uvw[,1],xlab='time [myr]',ylab='U [km/s]',type='l')
    plot(out[,1,1],uvw[,2],xlab='time [myr]',ylab='V [km/s]',type='l')
    plot(out[,1,1],uvw[,3],xlab='time [myr]',ylab='W [km/s]',type='l')
    plot(out[,1,1],sqrt(uvw[,1]^2+uvw[,2]^2+uvw[,3]^2),xlab='time [myr]',ylab='UVW [km/s]',type='l')
    dev.off()
    save(out,file=gsub('.pdf','.Robj',pdf.name))
}
