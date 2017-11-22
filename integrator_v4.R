#version final: adapt integrator.R to integrate TNO orbit
library("pracma")
library("numDeriv")
tide.HM<-function(t,state,par.tide){
with(as.list(c(state,par.tide)),{
  Omegar <- Omega-Omega0*t
  dL <- 0
  dG <- -5*L^2*(L^2-G^2)/(2*mu^2)*(cos(w)*sin(w)*(g3*(1-Theta^2/G^2)+(g1*sin(Omegar)^2+g2*cos(Omegar)^2)*Theta^2/G^2-g1*cos(Omegar)^2-g2*sin(Omegar)^2)-(g1-g2)*(cos(w)^2-sin(w)^2)*cos(Omegar)*sin(Omegar)*Theta/G)
  dTheta <- L^2*(g1-g2)/(2*mu^2)*(5*(L^2-G^2)*Theta/G*cos(w)*sin(w)*(cos(Omegar)^2-sin(Omegar)^2)+sin(Omegar)*cos(Omegar)*(G^2-Theta^2+5*(L^2-G^2)*(cos(w)^2-sin(w)^2*Theta^2/G^2)))
  dw <- L^2*G/(2*mu^2)*(g3*(1-5*sin(w)^2*(1-L^2*Theta^2/G^4))+(g1*cos(Omegar)^2+g2*sin(Omegar)^2)*(1-5*cos(w)^2)-5*(g1*sin(Omegar)^2+g2*cos(Omegar)^2)*L^2*Theta^2/G^4*sin(w)^2+5*(g1-g2)*cos(w)*sin(w)*cos(Omegar)*sin(Omegar)*(G^2+L^2)*Theta/G^3)
  dOmega <- L^2/(2*G*mu^2)*((g1*sin(Omegar)^2+g2*cos(Omegar)^2-g3)*(G^2+5*(L^2-G^2)*sin(w)^2)*Theta/G-5*(g1-g2)*(L^2-G^2)*cos(w)*sin(w)*cos(Omegar)*sin(Omegar))
  return(list(c(dw,dOmega,dL,dG,dTheta)))
})    
}

tide.CM<-function(t,state,par.tide){
with(as.list(c(state,par.tide)),{
  mu <- par.tide$mu
  g1 <- g1.func(t)
  g2 <- g2.func(t)
  g3 <- g3.func(t)
#  g1 <- par.tide$g1
#  g2 <- par.tide$g2
#  g3 <- par.tide$g3
  Omega0 <- par.tide$Omega0
  r <- sqrt(x^2+y^2+z^2)
  xp <- x*cos(Omega0*t)+y*sin(Omega0*t)
  yp <- -x*sin(Omega0*t)+y*cos(Omega0*t)
  dx <- xdot
  dy <- ydot
  dz <- zdot
  dxdot <- -mu*x/r^3-g1*xp*cos(Omega0*t)+g2*yp*sin(Omega0*t)
  dydot <- -mu*y/r^3-g1*xp*sin(Omega0*t)-g2*yp*cos(Omega0*t)
  dzdot <- -mu*z/r^3-g3*z
  return(list(c(dx,dy,dz,dxdot,dydot,dzdot)))
})    
}

map3 <- function(t,dt,state,par.tide){
  L <- state[1]
  G <- state[2]
  Theta <- state[3]
  w <- state[4]
  Omega <- state[5]
  Omega0 <- par.tide$Omega0
  mu <- par.tide$mu
  g1 <- par.tide$g1
  g2 <- par.tide$g2
  g3 <- par.tide$g3
#  g1 <- par.tide$g1
#  g2 <- par.tide$g2
#  g3 <- par.tide$g3
  Omegar <- Omega-Omega0*t
  dL <- 0
  dG <- -5*L^2*(L^2-G^2)/(2*mu^2)*(cos(w)*sin(w)*(g3*(1-Theta^2/G^2)+(g1*sin(Omegar)^2+g2*cos(Omegar)^2)*Theta^2/G^2-g1*cos(Omegar)^2-g2*sin(Omegar)^2)-(g1-g2)*(cos(w)^2-sin(w)^2)*cos(Omegar)*sin(Omegar)*Theta/G)
  dTheta <- L^2*(g1-g2)/(2*mu^2)*(5*(L^2-G^2)*Theta/G*cos(w)*sin(w)*(cos(Omegar)^2-sin(Omegar)^2)+sin(Omegar)*cos(Omegar)*(G^2-Theta^2+5*(L^2-G^2)*(cos(w)^2-sin(w)^2*Theta^2/G^2)))
  dw <- L^2*G/(2*mu^2)*(g3*(1-5*sin(w)^2*(1-L^2*Theta^2/G^4))+(g1*cos(Omegar)^2+g2*sin(Omegar)^2)*(1-5*cos(w)^2)-5*(g1*sin(Omegar)^2+g2*cos(Omegar)^2)*L^2*Theta^2/G^4*sin(w)^2+5*(g1-g2)*cos(w)*sin(w)*cos(Omegar)*sin(Omegar)*(G^2+L^2)*Theta/G^3)
  dOmega <- L^2/(2*G*mu^2)*((g1*sin(Omegar)^2+g2*cos(Omegar)^2-g3)*(G^2+5*(L^2-G^2)*sin(w)^2)*Theta/G-5*(g1-g2)*(L^2-G^2)*cos(w)*sin(w)*cos(Omegar)*sin(Omegar))
  dOmegar <- dOmega-Omega0
######calculate the second derivatives
  x <- c(L,G,Theta,w,Omegar)
  f1 <- function(x) c(dL, dG, dTheta, dw, dOmegar)
  f2 <- function(x) jacobian(f1,c(L,G,Theta,w, Omegar))%*%f1(x)
######calculate the third derivatives
  f3 <- function(x) jacobian(f2,c(L,G,Theta,w, Omegar))%*%f1(x)
  L <- L+dt*f1(x)[1]+dt^2/2*f2(x)[1]+dt^3/6*f3(x)[1]
  G <- G+dt*f1(x)[2]+dt^2/2*f2(x)[2]+dt^3/6*f3(x)[2]
  Theta <- Theta+dt*f1(x)[3]+dt^2/2*f2(x)[3]+dt^3/6*f3(x)[3]
  w <- w+dt*f1(x)[4]+dt^2/2*f2(x)[4]+dt^3/6*f3(x)[4]
  Omega <- Omega+dt*f1(x)[5]+dt^2/2*f2(x)[5]+dt^3/6*f3(x)[5]
  return(c(t+dt,w,Omega,L,G,Theta))
}

############################################################
##The following functions are for LARKS model of Oort cloud
############################################################
cv2ks <- function(x,y,z,vx,vy,vz,alpha){
  r <- sqrt(x^2+y^2+z^2)
  if(x>=0){
    u <- sqrt(alpha/(2*(r+x)))*c(0,r+x,y,z)
  }else{
    u <- sqrt(alpha/(2*(r-x)))*c(-z,y,r-x,0)}  
  U <- 2/alpha*c(u[1]*vx+u[4]*vy-u[3]*vz,u[2]*vx+u[3]*vy+u[4]*vz,-u[3]*vx+u[2]*vy-u[1]*vz,-u[4]*vx+u[1]*vy+u[2]*vz)
  val <- u[2]*U[1]-u[1]*U[2]-u[4]*U[3]+u[3]*U[4]
###give a test of this identity
  if(abs(val)>0.001){
     print("u1*U0-u0*U1-u3*U2+u2*U3!=0")}
  return(c(u=u,U=U))
}

ks2cv <- function(u,U,alpha){
  x <- (u[1]^2+u[2]^2-u[3]^2-u[4]^2)/alpha
  y <- 2*(u[2]*u[3]+u[1]*u[4])/alpha
  z <- 2*(u[2]*u[4]-u[1]*u[3])/alpha
  r <- sqrt(x^2+y^2+z^2)
  vx <- 1/(2*r)*(u[1]*U[1]+u[2]*U[2]-u[3]*U[3]-u[4]*U[4])
  vy <- 1/(2*r)*(u[4]*U[1]+u[3]*U[2]+u[2]*U[3]+u[1]*U[4])
  vz <- 1/(2*r)*(-u[3]*U[1]+u[4]*U[2]-u[1]*U[3]+u[2]*U[4])
  return(c(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz))
}

H0 <- function(x,y,z,vx,vy,vz,par.tide){
  mu <- par.tide$mu
  r <- sqrt(x^2+y^2+z^2)
  val <- 1/2*(vx^2+vy^2+vz^2)-mu/r
  return(val)
}

H1 <- function(t,x,y,z,par.tide){
  g1 <- par.tide$g1
  g2 <- par.tide$g2
  g3 <- par.tide$g3
  Omega0 <- par.tide$Omega0
  xp <- x*cos(Omega0*t)+y*sin(Omega0*t)
  yp <- -x*sin(Omega0*t) + y*cos(Omega0*t)
  val <- 1/2*(g1*xp^2+g2*yp^2+g3*z^2)
  return(val)
}

pH1.pu <- function(t,x,y,z,par.tide){
  r <- sqrt(x^2+y^2+z^2)
  if(x>=0){
    u <- sqrt(alpha/(2*(r+x)))*c(0,r+x,y,z)
  }else{
    u <- sqrt(alpha/(2*(r-x)))*c(-z,y,r-x,0)}  
  val <- grad(function(u) H1(t,x,y,z,par.tide),u)
}

K0 <- function(u,U,alpha,par.tide){
  x <- ks2cv(u,U,alpha)[1]
  y <- ks2cv(u,U,alpha)[2]
  z <- ks2cv(u,U,alpha)[3]
  vx <- ks2cv(u,U,alpha)[4]
  vy <- ks2cv(u,U,alpha)[5]
  vz <- ks2cv(u,U,alpha)[6]
  kappa0 <- H0(x,y,z,vx,vy,vz,par.tide)
  return(kappa0)
}

K1 <- function(u,t,U,alpha,par.tide){
  x <- ks2cv(u,U,alpha)["x"]
  y <- ks2cv(u,U,alpha)["y"]
  z <- ks2cv(u,U,alpha)["z"]
  kappa1 <- H1(t,x,y,z,par.tide)
  return(kappa1)
}

Us.ks <- function(t,u,U,alpha,par.tide){
  -K0(u,U,alpha,par.tide)-K1(u,t,U,alpha,par.tide)
}

Us.cv <- function(t,x,y,z,vx,vy,vz,par.tide){
  -H0(x,y,z,vx,vy,vz,par.tide)-H1(t,x,y,z,par.tide)
}

M0 <- function(u,U,alpha,par.tide){
  1/2*sum(U^2)+(4*Us/alpha^2)*sum(u^2)
}

M1 <- function(t,x,y,z,par.tide){
  4*sum(u^2)/alpha^2*H1(t,x,y,z,par.tide)
}

M1.ts <- function(t,u,U,alpha,par.tide){
  x <- ks2cv(u,U,alpha)[1]
  y <- ks2cv(u,U,alpha)[2]
  z <- ks2cv(u,U,alpha)[3]
  val <- M1(t,x,y,z,par.tide)
  return(val)
}

hess.M1 <- function(t,u,U,alpha,par.tide){
  x <- ks2cv(u,U,alpha)[1]
  y <- ks2cv(u,U,alpha)[2]
  z <- ks2cv(u,U,alpha)[3]
  val <- hessian(function(u) M1.ts(t,u,U,par.tide),u)
  return(val)
}
  
Phi0 <- function(u,t,U,Us,alpha,Delta,par.tide){
  if(Us>0){
    w <- 2*sqrt(2*Us)/alpha
    v <- u*cos(w*Delta)+U*w^-1*sin(w*Delta)
    V <- -u*w*sin(w*Delta)+U*cos(w*Delta)
    Vs <- Us
    t.new <- t+2*Delta/alpha^2*(sum(u^2)+sum(U^2)/w^2)+2*(u%*%U-v%*%V)/(alpha^2*w^2)
    C1 <- u%*%u + U%*%U/w^2
    C2 <- v%*%v + V%*%V/w^2
    err <- abs((C1-C2)/C1)
    if(err>0.01){
      print("The expression 'u^2+U^2*w^-2' is not invariant!")
      print(err)}
  }else{
      w <- 2*sqrt(-2*Us)/alpha
      v <- u*cosh(w*Delta)+U*w^-1*sinh(w*Delta)
      V <- u*w*sinh(w*Delta)+U*cosh(w*Delta)
      Vs <- Us
      t.new <- t+2*Delta/alpha^2*(sum(u^2)-sum(U^2)/w^2)-2*(u%*%U-v%*%V)/(alpha^2*w^2)
    C1 <- u%*%u - U%*%U/w^2
    C2 <- v%*%v - V%*%V/w^2
    err <- abs((C1-C2)/C1)
    if(err>0.01){
      print("The expression 'u^2+U^2*w^-2' is not invariant!")
      print(err)}
    }
  #GIVE A TEST FOR U AND u
  return(c(u=v,t=t.new,U=V,Us=Vs))
}

F <- function(u,t,U,alpha,par.tide){
  8*K1(u,t,U,alpha,par.tide)/alpha^2*u + 4*sum(u^2)/alpha^2*grad(function(u) K1(u,t,U,alpha,par.tide),u)
}

Fs <- function(u,t,U,alpha,par.tide){
  4*sum(u^2)/alpha^2*grad(function(t) K1(u,t,U,alpha,par.tide),t)
}

Phi1 <- function(u,t,U,Us,alpha,Delta,par.tide){
  v <- u
  t.new <- t
  V <- U-Delta*F(u,t,U,alpha,par.tide)
  Vs <- Us-Delta*Fs(u,t,U,alpha,par.tide)
  return(c(u=v,t=t.new,U=V,Us=Vs))
}

Phic <- function(u,t,U,Us,alpha,Delta,par.tide){
  v <- u
  t.new <- t
  V <- U-2*Delta*jacobian(function(u) F(u,t,U,alpha,par.tide), u)%*%F(u,t,U,alpha,par.tide)
  Vs <- Us-2*Delta*sum(jacobian(function(t) F(u,t,U,alpha,par.tide), t)*F(u,t,U,alpha,par.tide))
  return(c(u=v,t=t.new,U=V,Us=Vs))
}

Phih <- function(u,t,U,Us,alpha,h,par.tide){
  d1 <- h/12
  d2 <- 5/12*h
  c2 <- (1/2-sqrt(5)/10)*h
  c3 <- h/sqrt(5)
  q <- -h^3*(13-5*sqrt(5))/288
  D <- c(q,d1,c2,d2,c3,d2,c2,d1,q)
  funs <- list(Phic,Phi1,Phi0,Phi1,Phi0,Phi1,Phi0,Phi1,Phic)
  for(k in 1:3){
    Delta <- D[k]
    f <- funs[[k]]
    map <- as.numeric(f(u,t,U,Us,alpha,Delta,par.tide))
    u <- map[1:4]
    t <- map[5]
    U <- map[6:9]
    Us <- map[10]
  }
  return(c(u=u,t=t,U=U,Us=Us))
}

#Cartesian to Kepler
cs2oe <- function(cs){
  R <- cs[1:3]
  V <- cs[4:6]
  r <- sqrt(R%*%R)
  mu <- 6.67384e-11*1.98855e30*1e-9*149597871^-3*(3600*24)^2*365.242199^2
  v <- sqrt(V%*%V)
######equinoctial orbital element
  H <- cross(R,V)
  h <- sqrt(sum(H^2))
#  cat('dim(V)=',dim(V),'; dim(H)=',dim(H),'\n')
  E <- cross(V,H)/mu-R/r
  e <- sqrt(sum(E^2))#1 kepler: eccentricity
  p <- sum(h^2)/mu
  a <- h^2/(mu*(1-e^2))#2 kepler: semi-major axis
  i <- acos(H[3]/h)#3 kepler: inclination 0<=i<pi
  N <- cross(c(0,0,1),H)
  n <- sqrt(sum(N^2))
  if(N[2]<0){
    Omega <- 2*pi-acos(N[1]/n)
  }else{
    Omega <- acos(N[1]/n)
  }  #4 kepler: Right ascensioin of ascending node
  if(E[3]<0){
    w <- 2*pi-acos(sum(N*E)/(n*e))
  }else{
    w <- acos(sum(N*E)/(n*e))
  }#5 kepler: Argument of Perigee
  if(sum(R*V)<0){
    if(abs(E%*%R-e*r)<10^-3){
      theta <- 2*pi
    }else if(abs(E%*%R+e*r)<10^-3){
      theta <- pi}else{
    theta <- 2*pi-acos(sum(E*R)/(e*r))}
  }else{
    if(abs(E%*%R-e*r)<10^-3){
      theta <- 0
    }else if(abs(E%*%R+e*r)<10^-3){
      theta <- pi}else{
    theta <- acos(sum(E*R)/(e*r))}
  }
  if(theta<2*pi & theta > pi){
      Ea <- 2*pi-acos((e+cos(theta))/(1+e*cos(theta)))
  }else{
      Ea <- acos((e+cos(theta))/(1+e*cos(theta)))
  }# Eccentric anomaly
  M <- Ea-e*sin(Ea)#6 kepler: mean anomaly
  q <- a*(1-e)
  oe <- c(e,q,i,w,Omega,M)
  return(oe)
}
eps3 <- function(m,e,x){
    t1 <- cos(x)
    t2 <- -1+e*t1
    t3 <- sin(x)
    t4 <- e*t3
    t5 <- -x+t4+m
    t6 <- t5/(1/2*t5*t4/t2+t2)
    return(t5/((1/2*t3-1/6*t1*t6)*e*t6+t2))
}
KeplerStart3 <- function(m,e){
    t34 <- e^2
    t35 <- e*t34
    t33 <- cos(m)
    return(m+(-1/2*t35+e+(t34+3/2*t33*t35)*t33)*sin(m))
}
kep.murison2 <- function(m,e,tol=1e-6){
    Mnorm <- m%%(2*pi)
    E0 <- KeplerStart3(Mnorm,e)
    Ntry <- 1000
    for(k in 1:Ntry){
        E1 <- E0-eps3(Mnorm,e,E0)
        if(k==Ntry) 'Kepler solver failed to converge!\n'
        if(all(abs(E1-E0)<tol)) break()
        E0 <- E1
    }
    return(E1)
}

###hyperbolic Keplerian function
kep.hyper <- function(m,e,tol=1e-12){
    if(abs(m)>10){
        E0 <- asinh(m/e)
    }else{
        E0 <- m
    }
    Ntry <- 1000
    for(k in 1:Ntry){
        E1 <- E0+(m-e*sinh(E0)+E0)/(e*cosh(E0)-1)
        if(k==Ntry) cat('Kepler solver failed to converge!\n')
        if(all(abs(E1-E0)<tol)) break()
        E0 <- E1
    }
    return(E1)
}

#Kepler to Cartesian
oe2cs <- function(oe){
  #deriv Ea by solving Kepler's equation M=E-e*sinE for E
  e <- oe[1]
  q <- oe[2]
  i <- oe[3]
  w <- oe[4]
  Omega <- oe[5]
  M <- oe[6]
  a <- q/(1-e)
#  mu <- 6.67384e-11*1.98855e30*1e-9*149597871^-3*(3600*24)^2*365.242199^2#m^3/s^2
  mu <- 4*pi^2
  if(e<1){
      Ea <- kep.murison2(M,e)
  }else{
      Ea <- kep.hyper(M,e)
  }
  cat('True anomaly:',Ea*180/pi,'\n')
  P <- c(cos(w)*cos(Omega)-sin(w)*cos(i)*sin(Omega),cos(w)*sin(Omega)+sin(w)*cos(i)*cos(Omega),sin(w)*sin(i))
  Q <- c(-sin(w)*cos(Omega)-cos(w)*cos(i)*sin(Omega),-sin(w)*sin(Omega)+cos(w)*cos(i)*cos(Omega),sin(i)*cos(w))
  if(e<1){
      R <- a*(cos(Ea)-e)*P + a*sqrt(1-e^2)*sin(Ea)*Q#position
      Edot <- sqrt(mu/a^3)/(1-e*cos(Ea))
      V <- -a*sin(Ea)*Edot*P + a*sqrt(1-e^2)*cos(Ea)*Edot*Q#velocity
  }else{
      R <- a*(cosh(Ea)-e)*P - a*sqrt(e^2-1)*sinh(Ea)*Q#position
      Edot <- sqrt(-mu/a^3)/(1-e*cosh(Ea))
      V <- -a*sinh(Ea)*Edot*P + a*sqrt(e^2-1)*cosh(Ea)*Edot*Q#velocity
  }
  cs <- c(R,V)
  return(cs)
}

###The following function have something wrong which make ele2cs and cs2oe not inversible
ele2cs <- function(oe){
  e <- oe[1]
  q <- oe[2]
  inc <- oe[3]
  w <- oe[4]
  Omega <- oe[5]
  M <- oe[6]
  a <- q/(1-e)
  mu <- 6.67384e-11*1.98855e30*1e-9*149597871^-3*(3600*24)^2*365.242199^2
  Ea <- kep.murison2(M,e)
  nu <- 2*atan(((1+e)/(1-e))^0.5*tan(Ea/2))
  p <- a*(1-e^2)
  r <- a*(1-e*cos(Ea))
  h <- (mu*a*(1-e^2))^0.5
  X <- r*(cos(Omega)*cos(w+nu)-sin(Omega)*sin(w+nu)*cos(inc))
  Y <- r*(sin(Omega)*cos(w+nu)+cos(Omega)*sin(w+nu)*cos(inc))
  Z <- r*sin(inc)*sin(w+nu)
  Vx <- X*h*e*sin(nu)/r/p-h/r*(cos(Omega)*sin(w+nu)+sin(Omega)*cos(w+nu)*cos(inc))
  Vy <- Y*h*e/r/p*sin(nu)-h/r*(sin(Omega)*sin(w+nu)-cos(Omega)*cos(w+nu)*cos(inc))
  Vz <- Z*h*e/r/p*sin(nu+h/r*sin(inc)*cos(w+nu))
  cs <- c(X,Y,Z,Vx,Vy,Vz)
  return(cs)
}
#funciton for deciding the new step size of LARKS model
h1 <- function(a1,r,Nstep){
  a0 <- 50000
  h0 <- a0^(3/2)/Nstep
  if(a1>0){
  h1 <- h0*(a1/a0)^(3/2)}else{
    h1 <- h0*(a0^(3/2)/(abs(a1)^(1/2)*r))}
}

Uutest <- function(u,U){
  val <- u[2]*U[1]-u[1]*U[2]-u[4]*U[3]+u[3]*U[4]
  if(val==0){
    print("The values of u and U are correct.")
  }else{
    print("The values of u and U are wrong!")}
}

orbit.HM <- function(element,Nstep,par.tide){
  Nper <- length(element[,1])-1
  time.step <- rep(NA,Nstep*Nper+1)
  phase <- array(data=NA,dim=c(Nper*Nstep+1,7),dimnames=c("Times","phase"))
  per.acc <- element[,1]
  e <- element[,2]
  q <- element[,3]
  I <- element[,4]
  w <- element[,5]
  Omega <- element[,6]
  a <- q/(1-e)
  period <- diff(per.acc)
  phase[1,1] <- 0
  phase[1,2:7] <- k2c(a[1],e[1],I[1],Omega[1],w[1],M=0,par.tide)
  time.step[1] <- 0
  for(i in 1:Nper){
    t <- seq(per.acc[i],per.acc[i+1],by=period[i]/Nstep)
    t <- t[-1]
    time.step[((i-1)*Nstep+2):(i*Nstep+1)] <- t
    n <- 2*pi/period[i]
    M <- n*(t-per.acc[i])
    for(j in 1:Nstep){
      phase[(i-1)*Nstep+j+1,2:7] <- k2c(a[i],e[i],I[i],Omega[i],w[i],M=M[j],par.tide)
    }
  }
  phase[,1] <- time.step
  E <- H0(phase[,2],phase[,3],phase[,4],phase[,5],phase[,6],phase[,7],par.tide)+H1(phase[,1],phase[,2],phase[,3],phase[,4],par.tide)
  return(list("phase"=phase,"E"=E))
}

oe2de <- function(oe,par.tide){
  de <- rep(NA,5)
  e <- oe[1]
  q <- oe[2]
  I <- oe[3]
  w <- oe[4]
  Omega <- oe[5]
  a <- q/(1-e)
  L <- sqrt(par.tide$mu*a)
  G <- L*sqrt(1-e^2)
  Theta <- G*cos(I)
  de <- c(w, Omega, L, G, Theta)
  return(de)
}

de2oe <- function(de, par.tide){
  w <- de[1]
  Omega <- de[2]
  L <- de[3]
  G <- de[4]
  Theta <- de[5]
  oe <- rep(NA,5)
  e <- sqrt(1-G^2/L^2)
  a <- L^2/par.tide$mu
  q <- a*(1-e)
  I <- acos(Theta/G)
  M <- 0
  oe <- c(e,q,I,w,Omega,M)
  return(oe)
}

ac <- function(par.comet,e){
  inda <- par.comet$inda
  inde <- par.comet$inde
  return(10^inda*(1-e)^inde)
}

#######This function is used to calcuate the orbit of minro bodies in the solar system perturbed by the Galactic tide.
orbit.integrator <- function(OE,tp,par.comet,par.tide,hmmethod,csmethod){
  Nstep <- par.comet$Nstep
  age <- par.comet$age
  qc <- par.comet$qc
  rc <- par.comet$rc
  ql <- par.comet$ql
  OE <- OE
  tp <- tp
  Nper <- par.comet$Nper
#  cat('Nper=',Nper,'\n')
  element <- c(tp,OE)
  Ntar <- 0
  Nesc <- 0
  Ttar <- NA
  Tesc <- NA
  as <- OE[2]/(1-OE[1])#element at last step
  qs <- OE[2]
  rs <- 2*as-qs
  de <- as.numeric(oe2de(OE,par.tide))
  DE <- c(w=de[1], Omega=de[2],L=de[3], G=de[4], Theta=de[5])
  for(k in 1:Nper){
      if(tp>age) break()
###OE ==> DE
      tper <- c(tp,tp+as^(3/2))#new times
      tper <- c(tp,tp + as^(3/2))
###integrator
      out <- as.numeric(ode(func=tide.HM,y=DE,parms=par.tide,times=tper,method=hmmethod)[2,])  
      DE <- c(w=out[2], Omega=out[3],L=out[4], G=out[5], Theta=out[6])
###DE ==> OE
      OE <- de2oe(DE,par.tide)
      if(any(is.na(OE)) | any(is.na(DE))){
          cat('k=',k,'\n')
          cat('OE=',OE,'\n')
      }
      tp <- tper[2]
      element <- rbind(element,c(tp,OE))
  }
  return(list(element=element,Ttar=Ttar,Tesc=Tesc,Ntar=Ntar,Nesc=Nesc))
}

stellar2galactic <- function(dvs,l,b,a){
  bs <- -b
  if(l>=0 & l<=pi){
    ls <- l+pi
  }
  if(l>pi){
    ls <- l-pi
  }
  rotx <- matrix(data=c(1,0,0,0,cos(a),-sin(a),0,sin(a),cos(a)),nrow=3,ncol=3)
  roty <- matrix(data=c(cos(bs),0,sin(bs),0,1,0,-sin(bs),0,cos(bs)),nrow=3,ncol=3)
  rotz <- matrix(data=c(cos(ls),-sin(ls),0,sin(ls),cos(ls),0,0,0,1),nrow=3,ncol=3)
  dvg <- rotz%*%roty%*%rotx%*%dvs
  return(dvg)
}

galactic2stellar <- function(vg, ls, bs, a){
  rotx <- matrix(data=c(1,0,0,0,cos(a),sin(a),0,-sin(a),cos(a)),nrow=3,ncol=3)
  roty <- matrix(data=c(cos(bs),0,-sin(bs),0,1,0,sin(bs),0,cos(bs)),nrow=3,ncol=3)
  rotz <- matrix(data=c(cos(ls),sin(ls),0,-sin(ls),cos(ls),0,0,0,1),nrow=3,ncol=3)
  vs <- rotx%*%roty%*%rotz%*%vg
  return(vs)
}

####function to transfer orbital element in the ecliptic reference frame to Galactic coordinates
ele2bl <- function(oe){
   N <- as.integer(length(oe)/6)
   beta <- rep(NA,N)
   lambda <- rep(NA,N)
   for(i in 1:N){
       if(N==1){
           OE <- oe
       }else{
           OE <- oe[i,]
       }
      CS <- oe2cs(OE)
      bl <- xyz2bl(CS[1],CS[2],CS[3])#its ecliptic coord.: lambda, beta
      betalambda <- ecliptic2gal(b=bl[1],l=bl[2])
#      bl <- equatorial2gal(alpha=betalambda[2],delta=betalambda[1])
#      bl <- equ2gal(alpha=betalambda[2],delta=betalambda[1])
      beta[i] <- betalambda[1]
      lambda[i] <- betalambda[2]
   }
   return(cbind(beta,lambda))
}

##transform orbital elements from ecliptic reference frame to the galactic reference frame
ele2gal <- function(oe){
   N <- as.integer(length(oe)/6)
   oe.gal <- oe
   for(i in 1:N){
       if(N==1){
           OE <- oe
       }else{
           OE <- oe[i,]
       }
       CS <- oe2cs(OE)
       Rnorm <- sqrt(CS[1:3]%*%CS[1:3])
       Vnorm <- sqrt(CS[4:6]%*%CS[4:6])
#       cat('CS=',CS,'\n')
       betalambda.pos <- xyz2bl(CS[1],CS[2],CS[3])#its ecliptic coord.: lambda, beta
       betalambda.vel <- xyz2bl(CS[4],CS[5],CS[6])#
       bl.pos <- ecliptic2gal(lambda=betalambda.pos[2],beta=betalambda.pos[1])
       bl.vel <- ecliptic2gal(lambda=betalambda.vel[2],beta=betalambda.vel[1])
       Rvec <- as.numeric(bl2xyz(bl.pos[1],bl.pos[2]))
       Vvec <-  as.numeric(bl2xyz(bl.vel[1],bl.vel[2]))
       CS.gal <- c(Rnorm*Rvec,Vnorm*Vvec)
       if(N==1){
           oe.gal <- cs2oe(CS.gal)
       }else{
           oe.gal[i,] <- cs2oe(CS.gal)
       }
   }
   return(oe.gal)
}

###transform orbital elements from galactic reference frame to ecliptic reference frame
ele2ecl <- function(oe){
   N <- as.integer(length(oe)/6)
   oe.ecl <- oe
   for(i in 1:N){
       if(N==1){
           OE <- oe
       }else{
           OE <- oe[i,]
       }
       CS <- oe2cs(OE)
       Rnorm <- sqrt(CS[1:3]%*%CS[1:3])
       Vnorm <- sqrt(CS[4:6]%*%CS[4:6])
       bl.pos <- xyz2bl(CS[1],CS[2],CS[3])#its ecliptic coord.: lambda, beta
       bl.vel <- xyz2bl(CS[4],CS[5],CS[6])#
       bl.pos <- gal2ecliptic(b=bl.pos[1],l=bl.pos[2])
       bl.vel <- gal2ecliptic(b=bl.vel[1],l=bl.vel[2])
       Rvec <- Rnorm*as.numeric(bl2xyz(bl.pos[1],bl.pos[2]))
       Vvec <- Vnorm*as.numeric(bl2xyz(bl.vel[1],bl.vel[2]))
       CS.ecl <- c(Rvec,Vvec)
       if(N==1){
           oe.ecl <- cs2oe(CS.ecl)
       }else{
           oe.ecl[i,] <- cs2oe(CS.ecl)
       }
   }
   return(oe.ecl)
}
####cartesian coordinates to longitude and latitude in spherical coordinates
xyz2bl <- function(x,y,z){
  b <- rep(NA,length(x))
  l <- rep(NA,length(x))
  for(j in 1:length(x)){
    l.test <- atan(y[j]/x[j])
    if(x[j]<0){
      l[j] <- l.test+pi
    }else if((x[j]>0) & (y[j]<0)){
      l[j] <- l.test+2*pi
    }else{
      l[j] <- l.test
    }
    b[j] <- atan(z[j]/sqrt(x[j]^2+y[j]^2))
  }
  return(cbind(b,l))
}
bl2xyz <- function(b.rad,l.rad){
    x <- cos(b.rad)*cos(l.rad)
    y <- cos(b.rad)*sin(l.rad)
    z <- sin(b.rad)
    return(cbind(x,y,z))
}
###ref: http://physics.stackexchange.com/questions/88663/converting-between-galactic-and-ecliptic-coordinates
ecliptic2gal <- function(lambda,beta){
    rad2deg <- 180/pi
####epoch J2000
    lambdaG <- 180.01/rad2deg
    lambdaB <- 266.84/rad2deg
    betaG <- 29.80/rad2deg
    betaB <- -5.54/rad2deg
    BK <- 96.43/rad2deg
    z <- sin(betaG)*sin(beta)+cos(betaG)*cos(beta)*cos(lambda-lambdaG)
    y <- cos(beta)*sin(lambda-lambdaG)
    x <- cos(betaG)*sin(beta)-sin(betaG)*cos(beta)*cos(lambda-lambdaG)
    tmp <- xyz2bl(x,y,z)
    bt <- tmp[1]
    lt <- (BK-tmp[2])%%(2*pi)
###test
    test <- FALSE
    if(test){
        cat('lambda=',lambda*rad2deg,'; beta=',beta*rad2deg,'\n')
        cat('l=',lt*rad2deg,'; b=',bt*rad2deg,'\n')
        zp <- sin(betaG)*sin(bt)+cos(betaG)*cos(bt)*cos(BK-lt)
        yp <- cos(bt)*sin(BK-lt)
        xp <- cos(betaG)*sin(bt)-sin(betaG)*cos(bt)*cos(BK-lt)
        tmp <- xyz2bl(xp,yp,zp)
        betat <- tmp[1]
        lambdat <- tmp[2]+lambdaG
        cat('lambdat=',(lambdat*rad2deg)%%(360),'; betat=',betat*rad2deg,'\n')
    }
    return(c(bt,lt))
}


####galactic coordinates to ecliptic coordinates
gal2ecliptic <- function(b,l){
    rad2deg <- 180/pi
####epoch J2000
    lambdaG <- 180.01/rad2deg
    lambdaB <- 266.84/rad2deg
    betaG <- 29.80/rad2deg
    betaB <- -5.54/rad2deg
    BK <- 96.43/rad2deg
    zp <- sin(betaG)*sin(b)+cos(betaG)*cos(b)*cos(BK-l)
    yp <- cos(b)*sin(BK-l)
    xp <- cos(betaG)*sin(b)-sin(betaG)*cos(b)*cos(BK-l)
    tmp <- xyz2bl(xp,yp,zp)
    beta <- tmp[1]
    lambda <- tmp[2]+lambdaG
    return(cbind(beta,lambda))
}

#####galactic coordinates to equatorial coordinates
gal2equ <- function(b,l){
####epoch 2000
    rad2deg <- 180/pi
    alphaG <- 192.85/rad2deg
    alphaB <- 266.40/rad2deg
    deltaG <- 27.13/rad2deg
    deltaB <- -28.94/rad2deg
    BK <- 122.9/rad2deg
    zp <- sin(deltaG)*sin(b)+cos(deltaG)*cos(b)*cos(BK-l)
    yp <- cos(b)*sin(BK-l)
    xp <- cos(deltaG)*sin(b)-sin(deltaG)*cos(b)*cos(BK-l)
    tmp <- xyz2bl(xp,yp,zp)
    delta <- tmp[1]
    alpha <- (tmp[2]+alphaG)%%(2*pi)
    return(cbind(delta,alpha))
}

equ2gal <- function(alpha,delta){
    rad2deg <- 180/pi
    alphaG <- 192.85/rad2deg
    alphaB <- 266.40/rad2deg
    deltaG <- 27.13/rad2deg
    deltaB <- -28.94/rad2deg
    BK <- 122.9/rad2deg
    z <- sin(deltaG)*sin(delta)+cos(deltaG)*cos(delta)*cos(alpha-alphaG)
    y <- cos(delta)*sin(alpha-alphaG)
    x <- cos(deltaG)*sin(delta)-sin(deltaG)*cos(delta)*cos(alpha-alphaG)
    tmp <- xyz2bl(x,y,z)
    bt <- tmp[1]
    lt <- (BK-tmp[2])%%(2*pi)
    test <- FALSE
###test
    if(test){
        cat('alpha=',alpha*rad2deg,'; delta=',delta*rad2deg,'\n')
        cat('b=',bt*rad2deg,'; l=',lt*rad2deg,'\n')
        zp <- sin(deltaG)*sin(bt)+cos(deltaG)*cos(bt)*cos(BK-lt)
        yp <- cos(bt)*sin(BK-lt)
        xp <- cos(deltaG)*sin(bt)-sin(deltaG)*cos(bt)*cos(BK-lt)
        tmp <- xyz2bl(xp,yp,zp)
        deltat <- tmp[1]
        alphat <- (tmp[2]+alphaG)%%(2*pi)
        cat('alphat=',(alphat*rad2deg)%%360,'; deltat=',deltat*rad2deg,'\n')
    }
    return(c(bt,lt))
}



###ref: http://www2.astro.psu.edu/users/rbc/a501/c1_spherical_astronomy.pdf
equtorial2gal <- function(alpha,delta){
####epoch J1950
    alpha0 <- 282.25/rad2deg
    delta0 <- 62.6/rad2deg
    l0 <- 33/rad2deg
    x <- cos(delta)*cos(alpha-alpha0)
    y <- cos(delta)*sin(alpha-alpha0)*cos(delta0)+sin(delta)*sin(delta0)
    z <- sin(delta)*cos(delta0)-cos(delta)*sin(alpha-alpha0)*sin(delta0)
    tmp <- xyz2bl(x,y,z)
    bt <- tmp[1]
    lt <- (tmp[2]+l0)%%(2*pi)
    return(c(bt,lt))
}

#####
e2g.rec <- function(x,y,z){
    N <- length(x)
    r <- sqrt(x^2+y^2+z^2)
    rgal <- c()
    for(j in 1:N){
        da <- xyz2bl(x[j],y[j],z[j])
        bl <- equ2gal(da[2],da[1])
        rgal <- rbind(rgal,r[j]*bl2xyz(bl[1],bl[2]))
    }
    return(rgal)
}
#####velocity kick from stellar encounters using classical impulse approximation (Rickman 1976)
dv.gal <- function(b.target,b.sun,mstar,vph){
###only allow one target, but could be multiple encounters
###b.target, b.sun are in unit of au, mstar is in unit of Msun, vph is in unit of au/yr
    N <- round(length(b.target)/3)
    if(N==1){
        b.target <- matrix(b.target,nrow=1)
        b.sun <- matrix(b.sun,nrow=1)
    }
    bt2 <- b.target[,1]^2+b.target[,2]^2+b.target[,3]^2
    bs2 <- b.sun[,1]^2+b.sun[,2]^2+b.sun[,3]^2
    dvx <- 2*mu/vph*(b.target[,1]/bt2-b.sun[,1]/bs2) 
    dvy <- 2*mu/vph*(b.target[,2]/bt2-b.sun[,2]/bs2) 
    dvz <- 2*mu/vph*(b.target[,3]/bt2-b.sun[,3]/bs2) 
#    cat('length(dvx)=',length(dvx),'\n')
    return(cbind(dvx,dvy,dvz))
}
