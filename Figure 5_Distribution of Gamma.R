## This R script generates Figure 5.

library(MASS)

## Kappa values based on the Netherlands data
  Kf<-c(.159,.160,.205,.173,.207,.154,.121,.208,.191,.151,.133,.174,.144,.135,.166,.208)
  
  Km<-c(.221, .195, .221, .176, .179, .151, .179,.134,.181,.109,.116, .197,.139,.125,.145,.134)
  
  K<-c(Kf,Km)
  
## Kappa values based on the Australian data
  Af<-c(.0591,.2175,.2061,.0939,.122,.2351,.2043,.2232,.1281)
  Am<-c(.1778,.1782,.1836,.1554,.1251,.1508,.2187,.1612,.2049)
  Aall<-c(Af,Am)

## Combining the Netherlands and Australian data 
  KA=c(K,Aall)
  
  hist(KA, main="Histogram for Kappa", xlab="Kappa", border="black", probability=TRUE, breaks=10,xlim=c(0,0.3),  ylim=c(0,15))
  
## Fit a gamma distribution
  
  fit_G  <- fitdistr(KA, "gamma"); fit_G

## Getting the 2 parameters of the gamma distribution
  para1<-  fit_G$estimate[1]; para2<-  fit_G$estimate[2]
  
  xvals <- seq(0,0.3,0.001)

## Draw the density curve
  dens3 <- dgamma(xvals, para1,para2)
  
  lines(xvals, dens3,col="blue")
  
  pdf("rplot.pdf")
  hist(KA, main="Histogram for Kappa", xlab="Kappa", border="black", probability=TRUE, breaks=10,xlim=c(0,0.3),  ylim=c(0,15))
  lines(xvals, dens3,col="blue")
  dev.off()

    