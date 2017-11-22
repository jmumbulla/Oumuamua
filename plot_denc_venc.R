write.table(tab[,c('ra','dec')],file='pos.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
pdf('denc_venc.pdf',6,6)
size <- 1.2
par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=size,cex.axis=size,cex=size)
for(j in 1:3){
    if(j==1){
        index <- 1:nrow(val)
        col <- 'grey'
        plot(val[,'denc'],val[,'venc'],xlab=expression(d[enc]*" [pc]"),ylab=expression(v[enc]*" [km/s]"),log='y',pch=20,xlim=range(6,val[,'denc'],val[,'denc.5per']),ylim=range(val[,'venc'],val[,'venc.5per']),col=col)
    }else{
        if(j==2){
            col <- 'blue'
            HIP.la <- c(30920,37766,67092,80824,113020)
            ind <- match(HIP.la,val[,'HIP'])
            index <- ind[!is.na(ind)]
        }else if(j==3){            
            col <- 'red'
            hip.castor <- c(108706,85523)#Castor moving group
            hip.ic2391 <- c(25119)
            hip.hyades <- c(19255)#Hyades Moving Group
            hip.octans.near <- c(116132)#Octans-near association
            hip.ursa.minor <- c(57548)#Ursa Major Moving Group
            hip.multi <- c(82817, 19849, 71681, 71683, 70890)#GJ 644/643 System; ** STF 518ABC; alpha Cen A, B, C
            hip.group <- c(hip.castor,hip.ic2391,hip.hyades,hip.octans.near,hip.ursa.minor,hip.multi)
            ind <- match(hip.group,val[,'HIP'])
            index <- ind[!is.na(ind)]
        }else{
            col <- 'steelblue'
            index <- which(val[,'denc']<2 & val[,'venc']<15)
        }
        points(val[index,'denc'],val[index,'venc'],col=col,pch=20)
    }
    arrows(val[index,'denc.5per'],val[index,'venc'],val[index,'denc.95per'],val[index,'venc'],length=0.03,angle=90,code=3,col=col)
    arrows(val[index,'denc'],val[index,'venc.5per'],val[index,'denc'],val[index,'venc.95per'],length=0.03,angle=90,code=3,col=col)
}
#legend('topright',inset=c(-0.5,0),legend=c('Single and binary','Local association','Group and multiple','Single home candidates'),col=c('grey','blue','red','steelblue'),pch=20,xpd=NA)
#inset=c(-0.55,0)
legend('topright',legend=c('Single and binary','Local association','Group and multiple'),col=c('grey','blue','red'),pch=20,xpd=NA)
#HIP 104539, 17288, 103749
dev.off()
