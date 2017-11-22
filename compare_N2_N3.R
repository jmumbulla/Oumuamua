pdf('N2_N3_combined.pdf',8,8)
size <- 1.1
par(mfrow=c(2,2),mar=c(4.2,4.2,1,1),cex.lab=size,cex.axis=size,cex=size)
###T simulations
load('ind92_m13FALSE.Robj')
times <- seq(-5000,5000,by=0.01)
inds <- sort(sample(1:length(times),1000))
plot(times[inds],oes1[inds,1,1],xlab='Time [Myr]',ylab='e',type='l',ylim=c(0,1),xlim=c(-5000,5000))
legend('topright',legend='AB in T sim.',bty='n')

inds <- sort(sample(which(times>-300 & times< 800),1000))
plot(times[inds],oes2[inds,1,1],xlab='Time [Myr]',ylab='e',type='l',xlim=c(-300,800),ylim=c(0,1))
legend('topright',legend='ABC in T sim.',bty='n')

###TE simulations
times <- seq(-5000,5000,by=0.01)
inds <- sort(sample(1:length(times),1000))
load('ind92_m13.Robj')
plot(times[inds],oes1[inds,1,1],xlab='Time [Myr]',ylab='e',type='l',ylim=c(0,1),xlim=c(-5000,5000))
legend('topleft',legend='AB in TE sim.',bty='n')

times <- seq(-5000,5000,by=0.001)
inds <- sort(sample(which(times>-600 & times< 600),1000))
plot(times[inds],oes2[inds,1,1],xlab='Time [Myr]',ylab='e',type='l',xlim=c(-600,600),ylim=c(0,1))
legend('topright',legend='ABC in\nTE sim.',bty='n')
dev.off()
