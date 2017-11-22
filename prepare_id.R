library(MASS)
source('function_orbit_v6.R')
library(abind)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    ind11 <- as.integer(args[1])
}else{
    ind11 <- 0
}
ind22 <- ind11
dt <- 1e-2
tmin <- -1
tmax <- 0#Myr
dt0 <- -(2017-1600)/1e6#From 2016-Dec-12 to 1599-Dec-12
ind.max <- 2
sim.enc <- TRUE
LA <- TRUE#if only simulate perturbations from Galactic tide
#LA <- FALSE
sim.type <- 'tide'
#   sim.type <- 'point'
Nsamp <- 0#if Nsamp=0, only nominal orbits will be calculated
Norbit <- Nsamp+1#including the nominal orbit
ic.type <- 'rand'#rand, bd, ,ubd,tubd (random, bound, tide-unbound)
rate <- rate0 <- 80
if(ic.type!='rand'){
    ind.start <- 1
}else{
    ind.start <- 0
}
Dmax <- 2e-3#kpc, maximum perihelion
tol.bs <- 1e-9
method <- 'bs'#bulirsch-stoer algorithm
separate <- TRUE
if(!separate){
    Dt <- tmax-tmin
}else{
#    Dt <- min(1,tmax-tmin)
    Dt <- min(0.1,tmax-tmin)
}
if(sim.enc){
    Dt <- 1
}
barycenter <- FALSE
#ind.star <- c(0,0,154426, 152859, 154433)
#ind.star <- c(0,0,49486)
#ind.star <- c(0,0,154425 ,152858 ,154432   ,9450 , 49541)
#ind.star <- c(0,0,48394, 138446,  58414,  49541, 202909,   9450, 152858, 154425, 154432, 179627)
load('results/info_enc.Robj')
#ind11 <- 1
#ind22 <- 10
if(!LA){
#    ind.star <- c(0,0,out.all[ind11:ind22,'id'])
    if(all(ind11==0)){
        ind.star <- c(0,0)
    }else{
        ind.star <- c(0,0,ind11)
    }
}else{
    ind.star <- c(ind11)
}
#ind.star <- c(0,0)
N0 <- length(which(ind.star==0))
N1 <- length(which(ind.star>0))
#Ms <- c(1,0,41,23)
if(N1>0){
#    Ms <- c(1,0,rep(1,N1))
    Ms <- c(1,0,2.7)
}else{
    Ms <- c(1,0)
}
Nt <- length(ind.star)
if(LA) Ms <- rep(1,Nt)
if(ind11==0) Ms <- 1e-8
#e.g. 
#Rscript scattering_experiment.R 0.000001 0 0.1 1 TRUE 154432 152858 2.007 0.123
#Rscript scattering_experiment.R 0.01 0 100 0 FALSE 161508 161287 1 1
ind.bary <- 1:min(2,Nt)
set.seed(9000)
err <- 1
verbose <- FALSE
vv <- 0
test <- TRUE
Ncores <- 1
element <- TRUE
mas2deg <- 1e-3/3600
sun.static <- FALSE
####simulation set up
fix.sun <- TRUE
####load files
if(!LA){
#file <- 'data/combined7cat_v1_err1.csv'
file <- 'results/enc5pc.csv'
#if(!exists('tab')){
    tab <- read.csv(file,check.names=FALSE)
#}
}else{
#    tab <- read.csv('LA_data.csv',check.names=FALSE,sep=',')
    tab <- read.csv('results/enc5pc.csv',check.names=FALSE,sep=',')
#    tmin <- max(-100,min(3*tab[ind.star,'tenc'],-1))
}
###################################
####Part I: global parameters
###################################
rad2deg <- 180/pi
hr2rad <- 2*pi/24
kpcmyr2auyr <- 1e3*206265/1e6
kpc2au <- 1e3*206265
auyr2kms <- 4.74047
kpcmyr2kms <- kpcmyr2auyr*auyr2kms
#mu <- (2*pi)^2#G in unit of AU^3/yr^2/Msun
mu <-  6.67384e-11*1.98855e30*1e-9*149597871^-3*(3600*24)^2*365.242199^2
tspan <- c(tmin,tmax)
time.ap <- seq(0,tspan[2],by=dt)
time.bp <- seq(0,-tspan[1],by=dt)
if(all(time.ap==0)){
    time.ap <- c()
}
if(all(time.bp==0)){
    time.bp <- c()
}
if(length(time.bp)>0){
    times <- c(-rev(time.bp),time.ap)
}else{
    times <- time.ap
}
times <- unique(times)
######################################################################
####Part IV: Galactic parameters
######################################################################
###galactic parameters; see encounter_simulation.R for details
par.gal <- list(Md=7.9e10,Mb=1.4e10,Mh=6.98e11,bb=0.35,bd=0.25,bh=24.0,ad=3.55,A=0.014,B=-0.014,G=4.50e-12)#potential
parameters<-list(Md=7.9e10,Mb=1.4e10,Mh=6.98e11,bb=0.35,bd=0.25,bh=24.0,ad=3.55,A=0.014,B=-0.014,G=4.50e-12)
arm <- TRUE
if(arm){
    arm.par1 <- list(alpha=4.25,Rmin=3.48,thetamin=0.262,extent=6.0,dR=0.75,Vp=0.020,m=2,rho0=2.5e7,Rs=7,H=0.18)
    arm.par2 <- list(alpha=4.25,Rmin=3.48,thetamin=3.141+0.262,extent=6.0,dR=0.75,Vp=0.020,m=2,rho0=2.5e7,Rs=7,H=0.18)
    Rmax1<- arm.par1$Rmin+arm.par1$extent
    Rmin1<- arm.par1$Rmin
    Rmax2<- arm.par2$Rmin+arm.par2$extent
    Rmin2<- arm.par2$Rmin
}

######################################################################
####Part V: Solar parameters
######################################################################
###solar parameters
if(!sun.static){
    R0 <- 8.27#kpc
#    R0 <- 8.0#kpc
    Vc0 <- vc0 <- Vcirc(R0,par.gal)/kpcmyr2auyr#kpc/myr
    U0 <- 11.1/auyr2kms#au/yr
    V0 <- 12.24/auyr2kms#au/yr
    W0 <- 7.25/auyr2kms#au/yr
    Z0 <- 0.026#kpc
}else{
    R0 <- Vc0 <- U0 <- V0 <- W0 <- Z0 <- 0
}
vp.sun <- sqrt(U0^2+V0^2+W0^2)*auyr2kms
id.mean <- list(R=R0,Rdot=-U0/kpcmyr2auyr,phi=0,phidot=-(vc0+V0/kpcmyr2auyr)/R0,z=Z0,zdot=W0/kpcmyr2auyr)
ini<- c("R","Rdot","phi","phidot","z","zdot")#the 6D initial conditions
vary.ini<- factor(c(),ini)#used to find the varing initial conditions
sd<- list(R=0,Rdot=0,phi=0,phidot=0,z=0,zdot=0)
if(is.element("R",vary.ini)){sd$R<- 0.5}
if(is.element("Rdot",vary.ini)){sd$Rdot<- 0.00036}
if(is.element("phi",vary.ini)){sd$phi<- 0}
if(is.element("phidot",vary.ini)){sd$phidot<- 0.0030}
if(is.element("z",vary.ini)){sd$z<- 0.003}
if(is.element("zdot",vary.ini)){sd$zdot<- 0.00038}
Nvar<- length(vary.ini)
####
#state.sun.nominal <- c(R=id.mean$R,Rdot=id.mean$Rdot,phi=id.mean$phi,phidot=id.mean$phidot,z=id.mean$z,zdot=id.mean$zdot)
state.sun.nominal <- c(x1=id.mean$R,y1=0,z1=Z0,vx1=-U0/kpcmyr2auyr,vy1=-V0/kpcmyr2auyr-vc0,vz1=W0/kpcmyr2auyr)
if(Nsamp>0){
    state.sun<- cbind(x1=rnorm(Nsamp,id.mean$R,sd$R),y1=0,z=rnorm(Nsamp,Z0,sd$z),vx1=rnorm(Nsamp,-U0/kpcmyr2auyr,sd$Rdot),vy1=rnorm(Nsamp,-V0/kpcmyr2auyr,id.mean$R*sd$phidot)-vc0,vz1=rnorm(Nsamp,W0/kpcmyr2auyr,sd$zdot))
    state.sun <- rbind(state.sun.nominal,state.sun)
}else{
    state.sun <- t(as.data.frame(state.sun.nominal))
}
Rs <- R0#The distance from galactic center to the Sun
Zs <- Z0
rhot0 <- dent(Rs,Zs,parameters)
#Vc0 <- Vcirc(R0,par.gal)#in unit of au/yr
####
#Heliocentric Galactic coordinates of A/2017 U1 at epoch 1599-Dec-12:
# 1011.693 1982.131 684.5223 au
#phii <- state.sun[,4]*dt0
xi <- -1011.693/kpc2au+state.sun[,1]
yi <- -1982.13/kpc2au+state.sun[,2]
zi <- 684.5223/kpc2au+state.sun[,3]
#Columba/Carina young association/moving group; ref. https://arxiv.org/pdf/0808.3362.pdf; Torres et al. 2008
xe <- c(42,-14)*1e-3+state.sun[,1]
ye <- c(56,94)*1e-3+state.sun[,2]
ze <- c(-47,-33)*1e-3+state.sun[,3]
#Heliostatic velocity of A/2017 U1 in the Galactic frame at epoch 1599-Dec-12:
# -2.41048 -4.730554 -1.630141 au/yr or -11.42681 -22.42505 -7.727633 km/s
if(length(which(ind.star==0))>2){
    UVWi <- rbind(c(-2.41048, -4.730554, -1.630141),c(-13.2, -21.8,-5.9)/auyr2kms,c(-10.2, -23.0, -4.4)/auyr2kms)#au/yr
}else{
    UVWi <- t(matrix(c(-2.41048, -4.730554, -1.630141)))#au/yr
}
Ulsr <- UVWi[,1]+U0#au/yr
Vlsr <- UVWi[,2]+V0+Vc0*kpcmyr2auyr
Wlsr <- UVWi[,3]+W0
if(test){
    xi <- state.sun[,1]
    yi <- zi <- 0
    Wlsr <- Ulsr <- 0
    Vlsr <- Vc0*kpcmyr2auyr
}
vxi <- -Ulsr#au/yr
vyi <- -Vlsr 
vzi <- Wlsr
if(length(which(ind.star==0))>2){
statei <- cbind(x2=c(xi,xe),y2=c(yi,ye),z2=c(zi,ze),vx2=vxi/kpcmyr2auyr,vy2=vyi/kpcmyr2auyr,vz2=vzi/kpcmyr2auyr)
}else{
statei <- cbind(x2=xi,y2=yi,z2=zi,vx2=vxi/kpcmyr2auyr,vy2=vyi/kpcmyr2auyr,vz2=vzi/kpcmyr2auyr)
}
state.enc <- c()
if(N1>0){
    source('transform.R')
}
##total
if(!LA){
if(Nsamp==0){
    state <- t(matrix(c(as.numeric(state.sun),as.vector(t(statei)),state.enc)))
}else{
    state <- cbind(as.numeric(state.sun),as.vector(t(statei)),state.enc)
    rownames(state) <- NULL
}
}else{
    if(ind11==0){
        state <- statei
    }else{
        state <- state.enc
    }
}
