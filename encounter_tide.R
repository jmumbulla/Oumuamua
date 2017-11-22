Nstep <- ceiling(max(tref)/Dt)
yini <- state[k,]
for(j in 1:Nstep){
    inds <- which(tref>=(j-1)*Dt & tref<=(j*Dt))
    tt <- tref[inds]-(j-1)*Dt#times with respect to the initial time
    if(k==1){
####generate encounters
#####initial state
        x <- val[inds[1],1+(1:Nt)*6-5]
        y <- val[inds[1],1+(1:Nt)*6-4]
        z <- val[inds[1],1+(1:Nt)*6-3]
        vx <- val[inds[1],1+(1:Nt)*6-2]
        vy <- val[inds[1],1+(1:Nt)*6-1]
        vz <- val[inds[1],1+(1:Nt)*6]
        xb <- x[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
        yb <- y[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
        zb <- z[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
        vxb <- vx[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
        vyb <- vy[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
        vzb <- vz[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
        R.ref <- sqrt(xb^2+yb^2)
        z.ref <- zb
        r.ref <- c(xb,yb,zb)
        v.ref <- c(vxb,vyb,vzb)
        rhot <- dent(R.ref,z.ref,parameters)
        ratio <- rhot/rhot0[k]
        enc <- simEnc.fixFi(Dt,Dmax,r.ref,v.ref,rate=max(1,round(rate0*ratio)),FALSE)
        if(FALSE){
            cat('dim(enc)=',dim(enc),'\n')
            dev.new()
            hist(enc[,'cosk'])
            Sys.sleep(2)
            dev.off()
        }
        ii <- sort(enc[,1],index.return=TRUE)$ix
        enc <- enc[ii,]
        trel <- enc[,1]#relative time
        enc[,1] <- sign*(enc[,1]+(j-1)*Dt)#absolute time
        encs <- rbind(encs,enc)
    }else{
        enc <- encs[sign*encs[,1]>=(j-1)*Dt & sign*encs[,1]<=(j*Dt),]
        trel <- sign*(enc[,1])-(j-1)*Dt
    }
####generate a sequence of digitalized encounter times
    ind.seq <- 1
    ind.pair <- c()
    for(l in 1:length(trel)){
        ind <- which.min(abs(tt-trel[l]))
        ind.seq <- c(ind.seq,ind)
        ind.pair <- rbind(ind.pair,c(l,ind))
    }
    ind.seq <- c(ind.seq,length(tt))
    ind.uni <- unique(ind.seq)
    for(i in 1:(length(ind.uni)-1)){
        tseq <- tt[ind.uni[i]:ind.uni[i+1]]
        if(method!='bs'){
            out.ap <- ode(func=orbit,y=yini,parms=parameters,times=sign*(tseq-tseq[1]),method=method)
        }else{
            tmp <- bulirsch_stoer(f=orbit, t=sign*(tseq-tseq[1]), y0=matrix(yini,ncol=1),tol=tol.bs)
            out.ap <- cbind(tseq,tmp)
        }
        if(length(tseq)>2){
            pp <- cbind(sign*(tseq[-1]+(j-1)*Dt),out.ap[-1,2:ncol(out.ap)])
        }else{
            pp <- matrix(c(sign*(tseq[-1]+(j-1)*Dt),out.ap[-1,2:ncol(out.ap)]),nrow=1)
        }
        kk <- which(ind.pair[,2]==ind.uni[i+1])
        if(length(kk)>0){
            xx <- pp[nrow(pp),1+(1:Nt)*6-5]
            yy <- pp[nrow(pp),1+(1:Nt)*6-4]
            zz <- pp[nrow(pp),1+(1:Nt)*6-3]
            vxx <- pp[nrow(pp),1+(1:Nt)*6-2]
            vyy <- pp[nrow(pp),1+(1:Nt)*6-1]
            vzz <- pp[nrow(pp),1+(1:Nt)*6]
            tmp <- calc.kick(xs=xx,ys=yy,zs=zz,vxs=vxx,vys=vyy,vzs=vzz,Ms=Ms,ind.bary=ind.bary,enc=enc[kk,],verbose=verbose)
            if(verbose){
                cat('range(pp[nrow(pp),1+(1:Nt)*6-2:0])=',range(pp[nrow(pp),1+(1:Nt)*6-2:0])*kpcmyr2kms,'km/s\n')
            }
            pp[nrow(pp),1+(1:Nt)*6-2] <- pp[nrow(pp),1+(1:Nt)*6-2]+tmp[,1]#dvx
            pp[nrow(pp),1+(1:Nt)*6-1] <- pp[nrow(pp),1+(1:Nt)*6-1]+tmp[,2]#dvy
            pp[nrow(pp),1+(1:Nt)*6] <- pp[nrow(pp),1+(1:Nt)*6]+tmp[,3]#dvz
        }
        val[inds[(ind.uni[i]+1):ind.uni[i+1]],] <- pp
        yini <- pp[nrow(pp),-1]
    }
}
