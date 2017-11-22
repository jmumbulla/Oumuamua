###settings
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    nt <- as.integer(args[1])
    group <- as.logical(args[2])
}else{
    nt <- 1
    group <- TRUE
}
if(group){
   nts <- 10*(nt-1)+1:10
}
rad2deg <- 180/pi
hr2rad <- 2*pi/24
kpcmyr2auyr <- 1e3*206265/1e6
kpc2au <- 1e3*206265
auyr2kms <- 4.74047
kpcmyr2kms <- kpcmyr2auyr*auyr2kms
Nt <- 5000
###load the orbit of the Oumuamua
load('I1_orbitFALSE_tmin-100_tmax0_dt0.001.Robj')
t0 <- out[,1,1]
x0 <- out[,8,1]
y0 <- out[,9,1]
z0 <- out[,10,1]
vx0 <- out[,11,1]
vy0 <- out[,12,1]
vz0 <- out[,13,1]
###load the orbits of stars
###check which folders correspond to the time span for the target's orbit: 2196-2440
if(FALSE){
    Nt <- 5000
    dt <- 0.001
    time.bp <- seq(-100,0,by=dt)
    time.ap <- seq(0,100,by=dt)
    times <- c(time.bp,time.ap)
    ind1 <- which(times==-100)
    ind2 <- which(times==0)[1]
    Nsub <- ceiling(length(times)/Nt)
    ndiv1 <- ceiling(ind1/Nsub)
    ndiv2 <- ceiling(ind2/Nsub)
}

###load
Ndiv <- 1000
for(nt in nts){
   folder <- paste0('Nt',Nt,'nt',nt)
   val <- c()
for(j in 1:998){
    f1 <- paste0(folder,'/orbit240_1000div',j,'_Nsamp0_dt0.001_',Nt,'time',nt,'.Robj')
    nstar <- 240
    if(!file.exists(f1)){
        f1 <- paste0(folder,'/orbit208_1000div',j,'_Nsamp0_dt0.001_',Nt,'time',nt,'.Robj')
        nstar <- 208
    }
    load(f1)
    tref <- times[times<=0 & times>=-100]
    ind0 <- match(tref,t0)
    ind1 <- match(tref,times)
    id <- 240*(j-1)+1:nstar
    for(k in 1:length(tref)){
        i0 <- ind0[k]
        i1 <- ind1[k]
        d <- sqrt((star.out[i1,1,]-x0[i0])^2+(star.out[i1,2,]-y0[i0])^2+(star.out[i1,3,]-z0[i0])^2)
        ind <- which(d<5e-3)
        if(length(ind)>0){
            v <- sqrt((star.out[i1,4,ind]-vx0[i0])^2+(star.out[i1,5,ind]-vy0[i0])^2+(star.out[i1,6,ind]-vz0[i0])^2)
            tmp <- data.frame(tref[k],id[ind],d[ind]*1e3,v*kpcmyr2kms,star.out[i1,,ind,drop=FALSE])
            val <- rbind(val,tmp)
        }
    }
}
fobj <- paste0('Oumuamua/enc_nt',nt,'.Robj')
cat('output Robj:\n',fobj,'\n')
if(length(val)>0){
    if(is.null(dim(val))) val <- t(matrix(val))
    colnames(val) <- c('time.myr','id','d.pc','v.kms','x.kpc','y.kpc','z.kpc','vx.kpc.myr','vy.kpc.myr','vz.kpc.myr')
}
save(val,file=fobj)
}
