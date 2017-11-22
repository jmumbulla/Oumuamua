library(pracma)
library(splus2R)
#library(msProcess)
library(deSolve)
library(nls2)
library(stats)
####prepare initial conditions
source('prepare_id.R',local=TRUE)
indAB1 <- indC1 <- indAB <- indC <- rep(0,Norbit)
#####define some global variables for ode 
if(method!='bs'){
    if(sim.type=='combined'){
        orbit <- Orbit.combined
    }else if(sim.type=='point'){
        orbit <- Orbit.point
    }else if(sim.type=='tide'){
        orbit <- Orbit.tide
#        sim.enc <- FALSE
    }
}else{
    if(sim.type=='combined'){
        orbit <- Orbit2.combined
    }else if(sim.type=='point'){
        orbit <- Orbit2.point
    }else if(sim.type=='tide'){
        orbit <- Orbit2.tide
#        sim.enc <- FALSE
    }
}
parameters$simenc <- sim.enc
parameters$dt <- dt
parameters$enc <- c()
encs <- c()
out <- array(NA,dim=c(length(times),6*Nt+1,Norbit))
t1 <- proc.time()
for(k in 1:Norbit){
###set up
###simulation
    if(sim.enc){
        out[times==0,,k] <- c(0,state[k,])
        if(length(time.ap)>1){
            sign <- 1
            tref <- time.ap
            val <- array(NA,dim=c(length(tref),6*Nt+1))
            val[1,] <- c(0,state[k,])
            source('encounter_tide.R',local=TRUE)
            inds <- which(times>0)
            out[inds,,k] <- val[-1,]
        }
        if(length(time.bp)>1){
            sign <- -1
            tref <- time.bp
            val <- array(NA,dim=c(length(tref),6*Nt+1))
            val[1,] <- c(0,state[k,])
            source('encounter_tide.R',local=TRUE)
            inds <- which(times<0)
            val <- val[-1,]
            index <- sort(val[,1],index.return=TRUE)$ix
            out[inds,,k] <- val[index,]
        }
    }else{
        val <- c(0,state[k,])
        if(length(time.ap)>1){
            if(method!='bs'){
                out.ap <- ode(func=orbit,y=state[k,],parms=parameters,times=time.ap,method=method)
                val <- rbind(val,out.ap[-1,])
            }else{
                if(!separate){
                    tmp <- bulirsch_stoer(f=orbit, t=time.ap, y0=matrix(state[k,],ncol=1),tol=tol.bs)
                }else{
                    tmp <- bs.step(f=orbit, t=time.ap, y0=matrix(state[k,],ncol=1),Dt=Dt,tol=tol.bs)
                }
                val <- rbind(val,cbind(time.ap[-1],tmp[-1,]))
            }
        }
        if(length(time.bp)>1){
            if(method!='bs'){
                out.bp<- ode(func=orbit,y=state[k,],parms=parameters,times=-time.bp,method=method)
            }else{
                if(!separate){
                    tmp <- bulirsch_stoer(f=orbit, t=-time.bp, y0=matrix(state[k,],ncol=1),tol=tol.bs)
                }else{
                    tmp <- bs.step(f=orbit, t=-time.bp, y0=matrix(state[k,],ncol=1),Dt=Dt,tol=tol.bs)
                }
                out.bp <- cbind(-time.bp,tmp)
            }
            out.bp <- out.bp[-1,]
            inds <- nrow(out.bp):1
            out.bp <- out.bp[inds,]
            val <- rbind(out.bp,val)
        }
        out[,,k] <- val
    }
}
t2 <- proc.time()
dur <- (t2-t1)[3]/3600#hour
cat('duration:',format(dur,digit=3),'hour\n')
source('plot_scattering.R')
