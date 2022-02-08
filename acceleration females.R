rm(list = ls())
library(xtable)
library(MASS)

# based on uniroot.all function from package rootSolve
# the package is now archived
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



Bs=rep(0,9)
Cs=rep(0,9)
gams=rep(0,9)
thets=rep(0,9)
omegas=rep(0,9)
Ns=rep(0,9)
NSEs=rep(0,9)
x_Is=rep(0,9)
x_As=rep(0,9)
x_Bs=rep(0,9)
x_Cs=rep(0,9)
thetSEs=rep(0,9)


m=1958 #controls cohort, m=1959 equivalent to 1894 cohort
######HMD only########
####read in data####
setwd("D:/Study new/ACST4600/R_updated")
#setwd("D:/Study/ACST4600/Data&R")
nlddeath=read.csv("NLDdeath.txt",header=T)
nldpop=read.csv("NLDpop.txt",header=T)
continuous=read.csv("mortality_1986-2015_selection.csv",header=T)
cohort=read.table("fclt.txt",header=T)
attach(continuous)
attach(nlddeath)


for (m in 1958:1972) {
  
  #####smooth continuous####
  
  #setwd("D:/Study/ACST4600/Data&R")
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
  
  radix=100000
  dx=rep(0,50)
  for (i in 0:43) {
    exposures[i+1]=radix
    dx[i+1]=qx[i+1]*radix
    radix=radix-dx[i+1]
  }
  l108=radix
  #####choosing N####
  
  lps=rep(0,14)
  for (N in 85:98) {
    #####pareto part
    x=65:(N-1)
    x1=66:(N)
    x2=N:108
    x3=(N+1):109
    
    log.lik <- function(theta){
      alpha <- theta[1]
      delta <- theta[2]
      gam <- theta[3]
      out <- sum(dx[1:(N-65)]*(-(exp(alpha))/(exp(delta))*((exp(exp(delta)))^x-1)+log(1-exp(-(exp(alpha))/(exp(delta))*(exp(exp(delta)))^x*((exp(exp(delta)))-1)))))+exposures[N-64]*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^N-1))-l65*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^65-1))+sum(dx[(N-64):44]*log((1+gam*((x2-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)-(1+gam*((x3-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)))+l108*log((1+gam*((108-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam))
      return(out)
    }
    
    theta.start <- c(-11,-3,2) 
    out <- optim(theta.start, log.lik, hessian = FALSE, control = list(fnscale=-1), method="Nelder-Mead")
    beta.hat <- out$par
    beta.hat
    alpha=beta.hat[1]
    delta=beta.hat[2]
    B=exp(beta.hat[1])
    C=exp(exp(beta.hat[2]))
    gam=beta.hat[3]
    
    lps[N-84]=sum(dx[1:(N-65)]*(-(exp(alpha))/(exp(delta))*((exp(exp(delta)))^x-1)+log(1-exp(-(exp(alpha))/(exp(delta))*(exp(exp(delta)))^x*((exp(exp(delta)))-1)))))+exposures[N-64]*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^N-1))-l65*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^65-1))+sum(dx[(N-64):44]*log((1+gam*((x2-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)-(1+gam*((x3-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)))+l108*log((1+gam*((108-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam))
  }
  
  which.max(lps)
  N=85+which.max(lps)-1
  NSE=(1/((lps[which.max(lps)]-lps[which.max(lps)-1])-(lps[which.max(lps)+1]-lps[which.max(lps)])))^(1/2)
  
  
  ####refit with N=87####
  
  x=65:(N-1)
  x1=66:(N)
  x2=N:108
  x3=(N+1):109
  log.lik <- function(theta){
    alpha <- theta[1]
    delta <- theta[2]
    gam <- theta[3]
    out <- sum(dx[1:(N-65)]*(-(exp(alpha))/(exp(delta))*((exp(exp(delta)))^x-1)+log(1-exp(-(exp(alpha))/(exp(delta))*(exp(exp(delta)))^x*((exp(exp(delta)))-1)))))+exposures[N-64]*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^N-1))-l65*((-(exp(alpha))/(exp(delta)))*((exp(exp(delta)))^65-1))+sum(dx[(N-64):44]*log((1+gam*((x2-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)-(1+gam*((x3-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam)))+l108*log((1+gam*((108-N)/(1/((exp(exp(delta)))^N*exp(alpha)))))^(-1/gam))
    return(out)
  }
  
  theta.start <- c(-11,-3,2) 
  out <- optim(theta.start, log.lik, hessian = TRUE, control = list(fnscale=-1), method="Nelder-Mead")
  beta.hat <- out$par
  beta.hat
  B=exp(beta.hat[1])
  C=exp(exp(beta.hat[2]))
  gam=beta.hat[3]
  thet=1/(C^N*B)
  thetSE = c(-N*thet*log(C), -thet)%*%solve(-out$hessian[1:2,1:2])%*%c(-N*thet*log(C), -thet)
  
  omega=N-thet/gam
  
  Bs[m-1957]=B
  Cs[m-1957]=C
  gams[m-1957]=gam
  thets[m-1957]=thet
  omegas[m-1957]=omega
  Ns[m-1957]=N
  NSEs[m-1957]=NSE
  thetSEs[m-1957]=thetSE
  
  
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
  
  x_Is[m-1957]=x_I
  x_As[m-1957]=x_A$root
  x_Bs[m-1957]=x_B[2]
  x_Cs[m-1957]=x_C
  
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
  plot(y,hG1,main=paste0("Hazard function Slope Females ",m-65),xlab="x",ylab="h'(x)",col=1,type = 'l')
  lines(y,hP1, col=2)
  abline(v=N,lty=2)
  abline(v=x_A,lty=2)
  abline(v=maxL,lty=2)
  mtext("N",side = 1,at = N)
  mtext(expression(x[A]),side = 1,at = x_A$root)
  mtext(expression(N + theta/kappa),side = 1,at = maxL)
  legend(N,hG1[length(hG1)],c("Gompertz","GPD"), col=1:2,lty=1)
  
  
  plot(y1,hG,main=paste0("Hazard function Females ",m-65),xlab="x",ylab="h(x)",col=1,type = 'l')
  lines(y,hP, col=2)
  abline(v=N,lty=2)
  abline(v=x_A,lty=2)
  abline(v=x_B,lty=2)
  abline(v=maxL,lty=2)
  mtext("N",side = 1,at = N)
  mtext(expression(x[A]),side = 1,at = x_A$root)
  mtext(expression(x[B]),side = 1,at = x_B[2])
  mtext(expression(N + theta/kappa),side = 1,at = maxL)
  legend(65,hG[length(hG)],c("Gompertz","GPD"), col=1:2,lty=1)
  
  
  plot(y1,FG,main=paste0("CDF Females ",m-65),xlab="x",ylab="F(x)",col=1,type = 'l')
  lines(y,FP, col=2)
  abline(v=N,lty=2)
  abline(v=x_A,lty=2)
  abline(v=x_B,lty=2)
  abline(v=x_I,lty=2)
  abline(v=maxL,lty=2)
  mtext("N",side = 1,at = N)
  mtext(expression(x[A]),side = 1,at = x_A$root)
  mtext(expression(x[B]),side = 1,at = x_B[2])
  mtext(expression(x[I]),side = 1,at = x_I)
  mtext(expression(N + theta/kappa),side = 1,at = maxL)
  legend(65,1,c("Gompertz","GPD"), col=1:2,lty=1)
  
}


convert=cbind(Bs,Cs,thets,x_Is,abs(gams),Ns,x_As,x_Bs,x_Cs,omegas)
xtable(convert,digits=c(1,6,4,2,2,3,0,2,2,2,2))













