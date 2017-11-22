state.enc <- array(NA,dim=c(Norbit,6*(Nt-N0)))
xyzuvw <- array(NA,dim=c(Norbit,6*(Nt-N0)))
eobs.all <- obs.all <- c()
inds <- which(ind.star!=0)
for(kk in 1:length(inds)){
    j <- inds[kk]
    ind <- ind.star[j]
######################################################################
####Part II: coordiantes transformation for a sample of initial observables 
##Monte Carlo approach
######################################################################
#####default observables
    obs <- c(ra=tab[ind,'ra'],dec=tab[ind,'dec'],plx=tab[ind,'parallax'],pmra=tab[ind,'pmra'],pmdec=tab[ind,'pmdec'],rv=tab[ind,'rv'])
    eobs <- c(era=tab[ind,'ra_error']*mas2deg,edec=tab[ind,'dec_error']*mas2deg,eplx=tab[ind,'parallax_error'],epmra=tab[ind,'pmra_error'],epmdec=tab[ind,'pmdec_error'],erv=tab[ind,'erv'])
    err.basic <- as.numeric(c(tab[ind,c('ra_error','dec_error')]*mas2deg,tab[ind,c('parallax_error','pmra_error','pmdec_error','erv')]))
    cor.basic <- as.numeric(tab[ind,c('ra_dec_corr', 'ra_parallax_corr', 'ra_pmra_corr', 'ra_pmdec_corr', 'dec_parallax_corr', 'dec_pmra_corr', 'dec_pmdec_corr', 'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr')])
    cov.basic <- cor.basic*c(err.basic[1]*err.basic[2:5],err.basic[2]*err.basic[3:5],err.basic[3]*err.basic[4:5],err.basic[4]*err.basic[5])
    tmp <- matrix(c(0,cov.basic[1:4],0,
                    0,0,cov.basic[5:7],0,
                    0,0,0,cov.basic[8:9],0,
                    0,0,0,0,cov.basic[10],0,
                    0,0,0,0,0,0,
                    0,0,0,0,0,0
                    ),nrow=6,byrow=TRUE)
    cov.full <- tmp+t(tmp)+diag(err.basic^2)
    obs.all <- rbind(obs.all,obs)
    eobs.all <- rbind(eobs.all,eobs)
    ra.deg <- obs['ra']
    dec.deg <- obs['dec']
    plx.mas <- obs['plx']
    pm.ra.masyr <- obs['pmra']
    pm.dec.masyr <- obs['pmdec']
    rv.kms <- obs['rv']
    era <- eobs['ra']
    edec <- eobs['dec']
    eplx <- eobs['plx']
    epm.ra <- eobs['pmra']
    epm.dec <- eobs['pmdec']
    erv <- eobs['rv']
###sampling
    ras <- ra.deg
    decs <- dec.deg
    pm.ras <- pm.ra.masyr
    pm.decs <- pm.dec.masyr
    plxs <- plx.mas
    rvs <- rv.kms
###ref. https://arxiv.org/pdf/1504.06535.pdf
###Tokovinin et al. 2017
#    if(tab[ind,'HIP.ID']==51966) rvs <- 21.16
#    if(tab[ind,'HIP.ID']==113020) rvs <- -1.3
#    if(tab[ind,'HIP.ID']==104539) rvs <- -12
#    if(tab[ind,'HIP.ID']==107556) rvs <- -6.3
    if(Nsamp>0){
        val <- try(mvrnorm(2*Nsamp,mu=obs,Sigma=cov.full),TRUE)
        if(class(val)=='try-error') val <- mvrnorm(2*Nsamp,mu=obs,Sigma=diag(diag(cov.full)))
        ras <- c(ras,val[1:Nsamp,1])
        decs <- c(decs,val[1:Nsamp,2])
        plxs <- c(plxs,val[which(val[,3]>0),3][1:Nsamp])#plx is always positive
        pm.ras <- c(pm.ras,val[1:Nsamp,4])
        pm.decs <- c(pm.decs,val[1:Nsamp,5])
        rvs <- c(rvs,val[1:Nsamp,6])
    }
####coordinate transformation 
    bl <- equ2gal(ras/180*pi,decs/180*pi)
    if(is.null(dim(bl))){
        b <- bl[1]
        l <- bl[2]
    }else{
        b <- bl[,1]
        l <- bl[,2]
    }
    ds <- 1/plxs#kpc
    pos.gal<- hel2gal(ds,b,l,Rs,Zs)
    Xgc <- pos.gal$X
    Ygc <- pos.gal$Y
    Zgc<- pos.gal$Z
    Rgc <- pos.gal$R
    phigc <- pos.gal$phi
#transform velocity
#in helio-static frame
    vel <- Ve2hg(pm.ra=pm.ras,pm.de=pm.decs,rv=rvs,ra=ras/180*pi,de=decs/180*pi,d=ds)
    uvw <- vel$vhg
    vp <- vel$vp#au/yr
    vr <- rvs/auyr2kms#au/yr
    vtot <- sqrt(vp^2+vr^2)
    qest <- ds*1e3*vp/sqrt(vp^2+vr^2)#pc
    dest <- ds*1e3*vr/sqrt(vp^2+vr^2)#pc
    tenc <- -dest*206.265/sqrt(vp^2+vr^2)#kyr
#in Galactic-static frame
####to compare with M13
    Xgal <- -(Xgc-state.sun[,1])*1e3#pc
    Ygal <- -Ygc*1e3
    Zgal <- (Zgc-state.sun[,5])*1e3
    Us <- uvw[,1]*auyr2kms
    Vs <- uvw[,2]*auyr2kms
    Ws <- uvw[,3]*auyr2kms
###
    Ulsr <- uvw[,1]+U0#au/yr
    Vlsr <- uvw[,2]+V0
    Wlsr <- uvw[,3]+W0
    vx <- -Ulsr#au/yr
    vy <- -Vlsr - Vc0*kpcmyr2auyr
    vz <- Wlsr
    rgc0 <- rep(0,3)
    vgc0 <- rep(0,3)
    state.new <- cbind(cbind(Xgc,Ygc,Zgc),cbind(vx,vy,vz)/kpcmyr2auyr)
    xyzuvw.new <- cbind(cbind(Xgal,Ygal,Zgal),cbind(Us,Vs,Ws))
    state.enc[,6*kk-5:0] <- state.new
    xyzuvw[,6*kk-5:0] <- xyzuvw.new
    Rdot.gc <- vx*cos(phigc)+vy*sin(phigc)#au/yr
    phidot.gc <- (-vx*sin(phigc)+vy*cos(phigc))/Rgc#au/yr/kpc
    zdot.gc <- vz#kpc/Myr
    tmp <- c(R=Rgc,Rdot=Rdot.gc/kpcmyr2auyr,phi=phigc,phidot=phidot.gc/kpcmyr2auyr,z=Zgc,zdot=zdot.gc/kpcmyr2auyr)
#    state.star <- rbind(state.star,tmp)
#    colnames(state.star) <- names(id.mean)
}
