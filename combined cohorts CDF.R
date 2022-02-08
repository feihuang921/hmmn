devtools::install_github("u5838836/STLT/STLT")

library(STLT)
library(MASS)
rm(list=ls())

uniroot.all <- function(
  f,
  interval,
  lower = min(interval),
  upper = max(interval),
  tol = .Machine$double.eps^0.2,
  maxiter = 1000,
  trace = 0,
  n = 100,
  ...
) {
  xseq <- seq(lower, upper, len = n + 1)
  f_xseq <- f(xseq, ...)
  out <- xseq[which(f_xseq == 0)]
  f_sign <- f_xseq[1:n] * f_xseq[2:(n + 1)]
  i_range <- which(f_sign < 0)
  for (i in i_range)  {
    out <- c(
      out,
      uniroot(f, lower = xseq[i],
              upper = xseq[i + 1],
              maxiter = maxiter,
              tol = tol,
              trace = trace,
              ...)$root
    )
  }
  out
}



###########-------------- Females------------- ############

#setwd("D:/Study/ACST4600/R_updated")
nlddeath=read.csv("NLDdeath.txt",header=T)
nldpop=read.csv("NLDpop.txt",header=T)
continuous=read.csv("mortality_1986-2015_selection.csv",header=T)
cohort=read.table("fclt.txt",header=T)
attach(continuous)
attach(nlddeath)
qxs=NULL

for (m in 1958:1973) {
  
  #####smooth continuous####
  
  nlddeath=read.csv("NLDdeath.txt",header=T)
  nldpop=read.csv("NLDpop.txt",header=T)
  continuous=read.csv("mortality_1986-2015_selection.csv",header=T)
  cohort=read.table("fclt.txt",header=T)
  attach(continuous)
  attach(nlddeath)
  
  
  #females 92 in 1986 cohort
  dx=rep(0,50)
  for (i in 0:27) {
    dx[i+1]=sum(as.numeric(as.character(nlddeath[Age==(65+i)&Year==(m+i)&Sex=="f",12])))
  }
  for (i in 28:49) {
    dx[i+1]=length(continuous[ydie==(m+i)&floor(lifespan/365.25)==(i+65)&gender=="vrouwen",1])
  }
  
  #l100=nldpop[33348,12]
  l65=100000
  
  attach(nldpop)
  exposures=rep(0,50)
  for (i in 0:27) {
    exposures[i+1]=sum(as.numeric(as.character(nldpop[Age==(65+i)&Year==(m+i)&Sex=="f",12])))
  }
  for (i in 28:49) {
    exposures[i+1]=sum(dx[(i+1):50])
  }
  qx=dx/exposures
  attach(cohort)
  qx[1:28]=as.numeric(as.character(cohort[Year==(m-65),4][66:93]))
  detach(cohort)
  
  qxs=cbind(qxs,qx)
}

qx = apply(qxs,1,mean,na.rm=TRUE)
stltfit = stlt(65:112,qx,hessian = TRUE,radix = 100000)
coefs = stltfit$coefficients
B=coefs$B
C=coefs$C
gam=coefs$gamma
N=coefs$N
omega=stltfit$Omega
thet = 1/(B*C^N)
thetSE = c(-N*thet*log(C), -thet)%*%solve(-stltfit$Hess[1:2,1:2])%*%c(-N*thet*log(C), -thet)
thetSE

B1=B
C1=C
N1=N

kappa<-abs(gam)
kappa1=kappa
hfun <- function(x,kappa=kappa1, theta=thet, B=B1,C=C1,N=N1) {  
  log(B)+log(log(C))+x*log(C) -log(kappa)+2*log(theta) +2*log(1-kappa*(x-N)/theta) }

x_I= (log(log(C))-log(B))/log(C)  
x_I    ##  x_I is point of inflection of Gompertz
maxL<- N+thet/kappa 
maxL    ##  maxL is lifelength
x_A <- uniroot(function(x) hfun(x)-0.0,c(N-10, maxL-0.1))
x_A   ## x_A is acceleration age 
# par(las=1,bty="l")
# curve(hfun,from=85,to=maxL-0.1,xname="x")
# abline(v=x_A $root,lty=3)
hfun2 <- function(x,kappa=kappa1, theta=thet, B=B1,C=C1,N=N1) {  
  log(B)+x*log(C) +log(theta) +log(1-kappa*(x-N)/theta) }
x_B<- uniroot.all(function(x) hfun2(x)-0.0,c(N-10, maxL-0.1) )
x_B ## x_B is where hazards of Gomp and Pareto are equal (one value should be N).
#curve(hfun2,from=85,to=maxL-0.1,xname="x")
x_C<- N-1/log(C) +thet/kappa
x_C   ## x_C is where derivatives of log hazards are equal 

y = seq(N,omega-0.001,0.001)
y1 = seq(65,omega-0.001,0.001) 
hG1 = B*log(C)*C^y
hP1 = kappa/((thet-kappa*(y-N))^2)
hG = B*C^y1
hP = 1/(thet-kappa*(y-N))
FG = 1- exp(-B/log(C)*(C^y1-1))
FGN = exp(-B/log(C)*(C^N-1))
FP = 1 - FGN*(1-kappa*(y-N)/thet)^(1/kappa)

par(mfrow=c(1,3))

plot(y1,FG,main=NULL,xlab='',ylab="F(x)",col=1,type = 'l')
lines(y,FP, col=2)
abline(v=N,lty=2)
abline(v=x_A,lty=2)
abline(v=x_B,lty=2)
abline(v=x_I,lty=2)
abline(v=maxL,lty=2)
mtext("N",side = 1,at = N,cex = 0.7)
mtext(expression(x[A]),side = 1,at = x_A$root,cex = 0.7)
mtext(expression(x[B]),side = 1,at = x_B[2],cex = 0.7,padj = 0.4)
mtext(expression(x[I]),side = 1,at = x_I,cex = 0.7)
mtext(expression(U),side = 1,at = maxL,cex = 0.7)
legend(65,1,c("Gompertz","GPD"), col=1:2,lty=1)




plot(y1,hG,main=NULL,xlab='',ylab="h(x)",col=1,type = 'l')
lines(y,hP, col=2)
abline(v=N,lty=2)
abline(v=x_A,lty=2)
abline(v=x_B,lty=2)
abline(v=maxL,lty=2)
mtext("N",side = 1,at = N,cex = 0.7)
mtext(expression(x[A]),side = 1,at = x_A$root,cex = 0.7)
mtext(expression(x[B]),side = 1,at = x_B[2],cex = 0.7,padj = 0.4)
mtext(expression(U),side = 1,at = maxL,cex = 0.7)
legend(65,hG[length(hG)],c("Gompertz","GPD"), col=1:2,lty=1)



plot(y,hG1,main=NULL,xlab='',ylab="h'(x)",col=1,type = 'l')
lines(y,hP1, col=2)
abline(v=N,lty=2)
abline(v=x_A,lty=2)
abline(v=maxL,lty=2)
mtext("N",side = 1,at = N,cex = 0.7)
mtext(expression(x[A]),side = 1,at = x_A$root,cex = 0.7)
mtext(expression(U),side = 1,at = maxL,cex = 0.7)
legend(N,hG1[length(hG1)],c("Gompertz","GPD"), col=1:2,lty=1)

#CDF
par(mfrow=c(1,1))
plot(y1,FG,xlab='',ylab='',col=2,type = 'n')
lines(y1[1:33000],FG[1:33000],col=2)
lines(y,FP, col=3)
abline(v=N,lty=2)
abline(v=x_A,lty=2)
abline(v=x_B,lty=2)
abline(v=x_I,lty=2)
abline(v=maxL,lty=2)
mtext("N",side = 1,at = N)
mtext(expression(x[A]),side = 1,at = x_A$root)
mtext(expression(x[B]),side = 1,at = x_B[2],padj = 0.3)
mtext(expression(x[I]),side = 1,at = x_I)
mtext(expression(U),side = 1,at = maxL)
legend(65,1,c("KME", "Gompertz","GPD"), col=1:3,lty=1)


####-----KM-----####

FG65 = 1- exp(-B/log(C)*(C^65-1))
lines(65:111,c(FG65,1-(1-FG65)*cumprod((1-qx[1:46]))),lty=1,col=1)


#### Zoomed in ####

plot(y,FP,xlab='',ylab='',col=3,type = 'l')
lines(y1,FG, col=2)
abline(v=N,lty=2)
abline(v=x_A,lty=2)
abline(v=x_B,lty=2)
abline(v=x_I,lty=2)
abline(v=maxL,lty=2)
mtext("N",side = 1,at = N)
mtext(expression(x[A]),side = 1,at = x_A$root)
mtext(expression(x[B]),side = 1,at = x_B[2])
mtext(expression(x[I]),side = 1,at = x_I)
mtext(expression(U),side = 1,at = maxL)
legend(N,1,c("KME", "Gompertz","GPD"), col=1:3,lty=1,cex = 0.8)


####-----KM-----####

FG65 = 1- exp(-B/log(C)*(C^65-1))
lines(65:111,c(FG65,1-(1-FG65)*cumprod((1-qx[1:46]))),lty=1,col=1)




###########-------------- Males------------- ############

#setwd("D:/Study new/ACST4600/R_updated")
nlddeath=read.csv("NLDdeath.txt",header=T)
nldpop=read.csv("NLDpop.txt",header=T)
continuous=read.csv("mortality_1986-2015_selection.csv",header=T)
cohort=read.table("mclt.txt",header=T)
attach(continuous)
attach(nlddeath)
qxs=NULL

for (m in c(1958:1973)) { #1966
  
  #####smooth continuous####
  
  #setwd("D:/Study/ACST4600/Data&R")
  nlddeath=read.csv("NLDdeath.txt",header=T)
  nldpop=read.csv("NLDpop.txt",header=T)
  continuous=read.csv("mortality_1986-2015_selection.csv",header=T)
  cohort=read.table("mclt.txt",header=T)
  attach(continuous)
  attach(nlddeath)
  
  
  #females 92 in 1986 cohort
  dx=rep(0,50)
  for (i in 0:27) {
    dx[i+1]=sum(as.numeric(as.character(nlddeath[Age==(65+i)&Year==(m+i)&Sex=="m",12])))
  }
  for (i in 28:49) {
    dx[i+1]=length(continuous[ydie==(m+i)&floor(lifespan/365.25)==(i+65)&gender=="mannen",1])
  }
  detach(nlddeath)
  #l100=nldpop[33348,12]
  l65=50000
  
  attach(nldpop)
  exposures=rep(0,50)
  for (i in 0:27) {
    exposures[i+1]=sum(as.numeric(as.character(nldpop[Age==(65+i)&Year==(m+i)&Sex=="m",12])))
  }
  for (i in 28:49) {
    exposures[i+1]=sum(dx[(i+1):50])
  }
  detach(nldpop)
  qx=dx/exposures
  attach(cohort)
  qx[1:28]=as.numeric(as.character(cohort[Year==(m-65),4][66:93]))
  detach(cohort)
  
  qxs=cbind(qxs,qx)
}

qx = apply(qxs,1,mean,na.rm=TRUE)
stltfit = stlt(65:112,qx,endN = 98,censorAge = 105,hessian = TRUE,radix = 50000)
coefs = stltfit$coefficients
B=coefs$B
C=coefs$C
gam=coefs$gamma
N=coefs$N
omega=stltfit$Omega
thet = 1/(B*C^N)
thetSE = c(-N*thet*log(C), -thet)%*%solve(-stltfit$Hess[1:2,1:2])%*%c(-N*thet*log(C), -thet)
thetSE


B1=B
C1=C
N1=N

kappa<-abs(gam)
kappa1=kappa
hfun <- function(x,kappa=kappa1, theta=thet, B=B1,C=C1,N=N1) {  
  log(B)+log(log(C))+x*log(C) -log(kappa)+2*log(theta) +2*log(1-kappa*(x-N)/theta) }

x_I= (log(log(C))-log(B))/log(C)  
x_I    ##  x_I is point of inflection of Gompertz
maxL<- N+thet/kappa 
maxL    ##  maxL is lifelength
x_A <- uniroot(function(x) hfun(x)-0.0,c(N-10, maxL-0.1))
x_A   ## x_A is acceleration age 
# par(las=1,bty="l")
# curve(hfun,from=85,to=maxL-0.1,xname="x")
# abline(v=x_A $root,lty=3)
hfun2 <- function(x,kappa=kappa1, theta=thet, B=B1,C=C1,N=N1) {  
  log(B)+x*log(C) +log(theta) +log(1-kappa*(x-N)/theta) }
x_B<- uniroot.all(function(x) hfun2(x)-0.0,c(N-10, maxL-0.1) )
x_B ## x_B is where hazards of Gomp and Pareto are equal (one value should be N).
#curve(hfun2,from=85,to=maxL-0.1,xname="x")
x_C<- N-1/log(C) +thet/kappa
x_C   ## x_C is where derivatives of log hazards are equal 

y = seq(N,omega-0.01,0.01)
y1 = seq(65,omega-0.01,0.01) 
hG1 = B*log(C)*C^y
hP1 = kappa/((thet-kappa*(y-N))^2)
hG = B*C^y1
hP = 1/(thet-kappa*(y-N))
FG = 1- exp(-B/log(C)*(C^y1-1))
FGN = exp(-B/log(C)*(C^N-1))
FP = 1 - FGN*(1-kappa*(y-N)/thet)^(1/kappa)

par(mfrow=c(1,3))

plot(y1,FG,main=NULL,xlab='',ylab="F(x)",col=1,type = 'l')
lines(y,FP, col=2)
abline(v=N,lty=2)
abline(v=x_A,lty=2)
abline(v=x_B,lty=2)
abline(v=x_I,lty=2)
abline(v=maxL,lty=2)
mtext("N",side = 1,at = N,cex = 0.7)
mtext(expression(x[A]),side = 1,at = x_A$root,cex = 0.7)
mtext(expression(x[B]),side = 1,at = x_B[2],cex = 0.7)
mtext(expression(x[I]),side = 1,at = x_I,cex = 0.7)
mtext(expression(U),side = 1,at = maxL,cex = 0.7)
legend(65,1,c("Gompertz","GPD"), col=1:2,lty=1)




plot(y1,hG,main=NULL,xlab='',ylab="h(x)",col=1,type = 'l')
lines(y,hP, col=2)
abline(v=N,lty=2)
abline(v=x_A,lty=2)
abline(v=x_B,lty=2)
abline(v=maxL,lty=2)
mtext("N",side = 1,at = N,cex = 0.7)
mtext(expression(x[A]),side = 1,at = x_A$root,cex = 0.7)
mtext(expression(x[B]),side = 1,at = x_B[2],cex = 0.7)
mtext(expression(U),side = 1,at = maxL,cex = 0.7)
legend(65,hG[length(hG)],c("Gompertz","GPD"), col=1:2,lty=1)



plot(y,hG1,main=NULL,xlab='',ylab="h'(x)",col=1,type = 'l')
lines(y,hP1, col=2)
abline(v=N,lty=2)
abline(v=x_A,lty=2)
abline(v=maxL,lty=2)
mtext("N",side = 1,at = N,cex = 0.7)
mtext(expression(x[A]),side = 1,at = x_A$root,cex = 0.7)
mtext(expression(U),side = 1,at = maxL,cex = 0.7)
legend(N,hG1[length(hG1)],c("Gompertz","GPD"), col=1:2,lty=1)

#CDF
par(mfrow=c(1,1))
plot(y1,FG,xlab='',ylab='',col=2,type = 'n')
lines(y1[1:33000],FG[1:33000],col=2)
lines(y,FP, col=3)
abline(v=N,lty=2)
abline(v=x_A,lty=2)
abline(v=x_B,lty=2)
abline(v=x_I,lty=2)
abline(v=maxL,lty=2)
mtext("N",side = 1,at = N)
mtext(expression(x[A]),side = 1,at = x_A$root)
mtext(expression(x[B]),side = 1,at = x_B[2],padj = 0.3)
mtext(expression(x[I]),side = 1,at = x_I)
mtext(expression(U),side = 1,at = maxL)
legend(65,1,c("KME", "Gompertz","GPD"), col=1:3,lty=1)


####-----KM-----####

FG65 = 1- exp(-B/log(C)*(C^65-1))
lines(65:111,c(FG65,1-(1-FG65)*cumprod((1-qx[1:46]))),lty=1,col=1)


#### Zoomed in ####

plot(y,FP,xlab="",ylab='',col=3,type = 'l')
lines(y1,FG, col=2)
abline(v=N,lty=2)
abline(v=x_A,lty=2)
abline(v=x_B,lty=2)
abline(v=x_I,lty=2)
abline(v=maxL,lty=2)
mtext("N",side = 1,at = N)
mtext(expression(x[A]),side = 1,at = x_A$root)
mtext(expression(x[B]),side = 1,at = x_B[2])
mtext(expression(x[I]),side = 1,at = x_I)
mtext(expression(U),side = 1,at = maxL)
legend(N,1,c("KME", "Gompertz","GPD"), col=1:3,lty=1,cex = 0.8)


####-----KM-----####

FG65 = 1- exp(-B/log(C)*(C^65-1))
lines(65:111,c(FG65,1-(1-FG65)*cumprod((1-qx[1:46]))),lty=1,col=1)









