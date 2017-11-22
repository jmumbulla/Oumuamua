source('function_orbit_v6.R')
tab <- read.table('ephemeris.txt')
RA <- (tab[,3]+tab[,4]/60+tab[,5]/3600)*360/24#deg
DEC <- tab[,6]/abs(tab[,6])*(abs(tab[,6])+tab[,7]/60+tab[,8]/3600)#deg
ra <- RA[1]
dec <- DEC[1]
era <- sd(RA)
edec <- sd(DEC)
Nsamp <- 1e5
ras <- rnorm(Nsamp,ra,era)
decs <- rnorm(Nsamp,dec,edec)
bl <- equ2gal(ras/180*pi,decs/180*pi)
vg <- c(-11.42681, -22.42505, -7.727633)#km/s
rg <- c(1011.693, 1982.131, 684.5223)#au
Vinf <- as.numeric(sqrt(vg%*%vg))
R <- as.numeric(sqrt(rg%*%rg))
uvw <- -Vinf*bl2xyz(bl[,1],bl[,2])
xyz <- R*bl2xyz(bl[,1],bl[,2])
cat('u',uvw[1,1],'\n')
cat('v',uvw[1,2],'\n')
cat('w',uvw[1,3],'\n')
cat('sd(u)',sd(uvw[,1]),'\n')
cat('sd(v)',sd(uvw[,2]),'\n')
cat('sd(w)',sd(uvw[,3]),'\n')
cat('Vinf=',mean(Vinf),'\n\n')
cat('x',xyz[1,1],'\n')
cat('y',xyz[1,2],'\n')
cat('z',xyz[1,3],'\n')
cat('sd(x)',sd(xyz[,1]),'\n')
cat('sd(y)',sd(xyz[,2]),'\n')
cat('sd(z)',sd(xyz[,3]),'\n')
ev <- sqrt(sd(uvw[,1])^2+sd(uvw[,2])^2+sd(uvw[,3])^2)
cat('sd(Vinf)=',ev,'km/s\n')
