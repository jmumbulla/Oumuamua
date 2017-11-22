Mt <- sum(Ms[ind.bary])
Mmu <- prod(Ms[ind.bary])/Mt
x <- state[,(0:(Nt-1))*6+1,drop=FALSE]
y <- state[,(0:(Nt-1))*6+2,drop=FALSE]
z <- state[,(0:(Nt-1))*6+3,drop=FALSE]
vx <- state[,(0:(Nt-1))*6+4,drop=FALSE]
vy <- state[,(0:(Nt-1))*6+5,drop=FALSE]
vz <- state[,(0:(Nt-1))*6+6,drop=FALSE]
###
xb <- x[,ind.bary,drop=FALSE]%*%Ms[ind.bary]/Mt#barycenter
yb <- y[,ind.bary]%*%Ms[ind.bary]/Mt
zb <- z[,ind.bary]%*%Ms[ind.bary]/Mt
vxb <- vx[,ind.bary]%*%Ms[ind.bary]/Mt#barycentric velocity
vyb <- vy[,ind.bary]%*%Ms[ind.bary]/Mt
vzb <- vz[,ind.bary]%*%Ms[ind.bary]/Mt
###First test whether B is bounded to A
dx <- (x-c(xb))*kpc2au
dy <- (y-c(yb))*kpc2au
dz <- (z-c(zb))*kpc2au
dvx <- (vx-c(vxb))*kpcmyr2auyr
dvy <- (vy-c(vyb))*kpcmyr2auyr
dvz <- (vz-c(vzb))*kpcmyr2auyr
r <- sqrt(dx^2+dy^2+dz^2)#au
dr <- sqrt((dx-dx[,1])^2+(dy-dy[,1])^2+(dz-dz[,1])^2)#au
v <- sqrt(dvx^2+dvy^2+dvz^2)#au/yr
dv <- sqrt((dvx-dvx[,1])^2+(dvy-dvy[,1])^2+(dvz-dvz[,1])^2)#au/yr
#dv <- 0.1296148/auyr2kms
##
bary <- FALSE
if(bary){
    V <- -(2*pi)^2*Mmu*Mt/dr
    K <- 0.5*Mmu*dv^2
    vesc <- sqrt(-V/(0.5*Mmu))
}else{
    V <- rep(0,Norbit)
    for(j in ind.bary){
        R <- sqrt((dx[,j]-dx[,-j])^2+(dy[,j]-dy[,-j])^2+(dz[,j]-dz[,-j])^2)
        if(is.null(dim(R))) R <- matrix(R,ncol=1)
        V <- V-(2*pi)^2*rowSums(Ms[j]*Ms[-j]/R)
    }
    K <- 0.5*v[,ind.bary]^2%*%Ms[ind.bary]
    vesc <- c()
    for(j in ind.bary){
        v2 <- (-V[1]-0.5*sum(v[1,-j]^2*Ms[-j]))/(0.5*Ms[j])
        if(v2>0){
            vesc <- c(vesc,sqrt(v2))
        }else{
            vesc <- c(vesc,0)
        }
    }
}
Energy <- K+V
cat('Energy=',Energy[1],'\n')
inds <- which(Energy<0)
indAB[inds] <- 1
cat('Ratio of bounded orbits of B to A:',length(inds)/Norbit,'\n')
cat('Ratio of unbounded orbits of B to A:',1-length(inds)/Norbit,'\n')
#source('useful_quantity.R')
#####binary
Mt <- sum(Ms[ind.bary])
xb <- x[,ind.bary]%*%Ms[ind.bary]/Mt
yb <- y[,ind.bary]%*%Ms[ind.bary]/Mt
zb <- z[,ind.bary]%*%Ms[ind.bary]/Mt
vxb <- vx[,ind.bary]%*%Ms[ind.bary]/Mt
vyb <- vy[,ind.bary]%*%Ms[ind.bary]/Mt
vzb <- vz[,ind.bary]%*%Ms[ind.bary]/Mt
###
if(Nt>2){
    ind <- 1
    dx1 <- (x[,3]-x[ind,1])*kpc2au
    dy1 <- (y[,3]-y[ind,1])*kpc2au
    dz1 <- (z[,3]-z[ind,1])*kpc2au
    dvx1 <- (vx[,3]-vx[ind,1])*kpcmyr2auyr
    dvy1 <- (vy[,3]-vy[ind,1])*kpcmyr2auyr
    dvz1 <- (vz[,3]-vz[ind,1])*kpcmyr2auyr
    dx2 <- (x[,3]-x[ind,2])*kpc2au
    dy2 <- (y[,3]-y[ind,2])*kpc2au
    dz2 <- (z[,3]-z[ind,2])*kpc2au
    dvx2 <- (vx[,3]-vx[ind,2])*kpcmyr2auyr
    dvy2 <- (vy[,3]-vy[ind,2])*kpcmyr2auyr
    dvz2 <- (vz[,3]-vz[ind,2])*kpcmyr2auyr
###
    dxb <- (x[,3]-c(xb))*kpc2au
    dyb <- (y[,3]-c(yb))*kpc2au
    dzb <- (z[,3]-c(zb))*kpc2au
    dvxb <- (vx[,3]-c(vxb))*kpcmyr2auyr
    dvyb <- (vy[,3]-c(vyb))*kpcmyr2auyr
    dvzb <- (vz[,3]-c(vzb))*kpcmyr2auyr
###
    d1 <- sqrt(dx1^2+dy1^2+dz1^2)
    v1 <- sqrt(dvx1^2+dvy1^2+dvz1^2)
    d2 <- sqrt(dx2^2+dy2^2+dz2^2)
    v2 <- sqrt(dvx2^2+dvy2^2+dvz2^2)
    V1 <- -(2*pi)^2*Ms[1]/d1
    V2 <- -(2*pi)^2*Ms[2]/d2
    V <- V1+V2
    K <- 0.5*Ms[1]*(dvxb^2+dvyb^2+dvzb^2)
    Energy <- K+V
    inds <- which(Energy<0)
    indC[inds] <- 1
###
    cat('Ratio of bounded orbits of C to A+B:',length(inds)/Norbit,'\n')
    cat('Ratio of unbounded orbits of C to A+B:',1-length(inds)/Norbit,'\n')
}
