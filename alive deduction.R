# setwd("D:/Study/ACST4600/Data&R")
# nlddeath=read.csv("NLDdeath.txt",header=T)
# nldpop=read.csv("NLDpop.txt",header=T)
# continuous=read.csv("mortality_1986-2015_selection.csv",header=T)
# attach(continuous)
# deathtimes=continuous[ydie==2015&gender=="mannen",7]/365.25
# library(xtable)
# 
# #female
# attach(nldpop)
# alivef=matrix(nrow = 30,ncol = 8)
# rownames(alivef)=1986:2015
# colnames(alivef)=92:99
# for (i in 1986:2015) {
#   for (j in 92:99) {
#     alivef[i-1985,j-91]=nldpop[Age==j&Year==i&Sex=="f",12]
#   }
# }
# 
# 
# #male
# attach(nldpop)
# alivem=matrix(nrow = 30,ncol = 8)
# rownames(alivem)=1986:2015
# colnames(alivem)=92:99
# for (i in 1986:2015) {
#   for (j in 92:99) {
#     alivem[i-1985,j-91]=nldpop[Age==j&Year==i&Sex=="m",12]
#   }
# }
# 
# xtable(alivem)
# xtable(alivef)

library(xtable)

############### Female ###############
####read in data####
m=1958 #controls cohort, m=1959 equivalent to 1894 cohort
nlddeath=read.csv("NLDdeath.txt",header=T)
nldpop=read.csv("NLDpop.txt",header=T)
continuous=read.csv("mortality_1986-2015_selection.csv",header=T)
cohort=read.table("fclt.txt",header=T)
attach(continuous)
attach(nlddeath)

Ns=c(97,98,98,98,97,98,98,98,97,98,98,96,97,98,97,95)
alive_at_N=rep(0,16)
alive_at_98=rep(0,16)
alive_at_100=rep(0,16)
alive_at_105=rep(0,16)
alive_at_110=rep(0,16)
alive_at_115=rep(0,16)

for (m in 1957:1972) {
  
  #females 92 in 1986 cohort
  dx=rep(0,50)
  for (i in 0:27) {
    dx[i+1]=sum(as.numeric(as.character(nlddeath[Age==(65+i)&Year==(m+i)&Sex=="f",12])))
  }
  for (i in 28:49) {
    dx[i+1]=length(continuous[ydie==(m+i)&floor(lifespan/365.25)==(i+65)&gender=="vrouwen",1])
  }
  
  #l100=nldpop[33348,12]
  #l65=1000000
  
  attach(nldpop)
  exposures=rep(0,50)
  for (i in 0:27) {
    exposures[i+1]=sum(as.numeric(as.character(nldpop[Age==(65+i)&Year==(m+i)&Sex=="f",12])))
  }
  for (i in 28:49) {
    exposures[i+1]=sum(dx[(i+1):50])
  }
  
  alive_at_N[m-1956] = exposures[Ns[m-1956]-64]
  alive_at_98[m-1956] = exposures[34]
  alive_at_100[m-1956] = exposures[36]
  alive_at_105[m-1956] = exposures[41]
  alive_at_110[m-1956] = exposures[46]
  alive_at_115[m-1956] = exposures[50]
}

alive_at_N



############### Male ###############
####read in data####
m=1958 #controls cohort, m=1959 equivalent to 1894 cohort
setwd("D:/Study/ACST4600/Data&R")
nlddeath=read.csv("NLDdeath.txt",header=T)
nldpop=read.csv("NLDpop.txt",header=T)
continuous=read.csv("mortality_1986-2015_selection.csv",header=T)
cohort=read.table("fmclt.txt",header=T)
attach(continuous)
attach(nlddeath)

Ns=c(89,91,89,98,98,98,95,96,94,97,98,93,98,95,98,96)
alive_at_N=rep(0,16)
alive_at_97=rep(0,16)
alive_at_100=rep(0,16)
alive_at_105=rep(0,16)
alive_at_110=rep(0,16)
alive_at_115=rep(0,16)

for (m in 1957:1972) {
  
  #females 92 in 1986 cohort
  dx=rep(0,50)
  for (i in 0:27) {
    dx[i+1]=sum(as.numeric(as.character(nlddeath[Age==(65+i)&Year==(m+i)&Sex=="m",12])))
  }
  for (i in 28:49) {
    dx[i+1]=length(continuous[ydie==(m+i)&floor(lifespan/365.25)==(i+65)&gender=="mannen",1])
  }
  
  #l100=nldpop[33348,12]
  #l65=1000000
  
  attach(nldpop)
  exposures=rep(0,50)
  for (i in 0:27) {
    exposures[i+1]=sum(as.numeric(as.character(nldpop[Age==(65+i)&Year==(m+i)&Sex=="m",12])))
  }
  for (i in 28:49) {
    exposures[i+1]=sum(dx[(i+1):50])
  }
  
  alive_at_N[m-1956] = exposures[Ns[m-1956]-64]
  alive_at_97[m-1956] = exposures[33]
  alive_at_100[m-1956] = exposures[36]
  alive_at_105[m-1956] = exposures[41]
  alive_at_110[m-1956] = exposures[46]
  alive_at_115[m-1956] = exposures[50]
}

alive_at_N

