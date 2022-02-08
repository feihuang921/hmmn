library(readr)
library(tidyverse)
library(lubridate)
library(DescTools)
#Import density data
densities=read_delim("areas_f_new_20210910.txt",col_names=F,delim=',') 

area_110=densities$X17[3]
area_115=densities$X17[4]
area_120=densities$X17[5]
area_125=densities$X17[6]

prob_115=area_115/area_110
prob_120=area_120/area_110
prob_125=area_125/area_110



## Import IDL data from "International Data Base on Longevity"
idl_complete <- read_delim("idl_complete.csv", 
                           ";", escape_double = FALSE, trim_ws = TRUE)

idl_complete$date_of_death=dmy(sub("/","-",sub("/","-",idl_complete$DDATE)))
idl_complete$date_of_birth=dmy(sub("/","-",sub("/","-",idl_complete$BDATE)))
idl_complete$YOB=year(idl_complete$date_of_birth)
idl_complete$YOD=year(idl_complete$date_of_death)
idl_complete[is.na(idl_complete$date_of_death),]$YOD=as.integer(idl_complete[is.na(idl_complete$date_of_death),]$DDATE)
idl_complete[is.na(idl_complete$date_of_birth),]$YOB=as.integer(idl_complete[is.na(idl_complete$date_of_death),]$BDATE)

nsupercents_by_year_and_country=idl_complete %>% filter(VALIDATION=="YES") %>% 
  filter(AGEYEARS>=110) %>% 
  group_by(DCOUNTRY,YOD) %>%
  summarize(n=n())


#Let's try modelling growth in all countries with >=10 supercentenarians
#BEL  11
#ESP  23
#EW   40
#FRA  28
#JPN  10
#USA  31

BEL_IDL=nsupercents_by_year_and_country[nsupercents_by_year_and_country$DCOUNTRY=="BEL",]
ESP_IDL=nsupercents_by_year_and_country[nsupercents_by_year_and_country$DCOUNTRY=="ESP",]
EW_IDL=nsupercents_by_year_and_country[nsupercents_by_year_and_country$DCOUNTRY=="EW",]
FRA_IDL=nsupercents_by_year_and_country[nsupercents_by_year_and_country$DCOUNTRY=="FRA",]
JPN_IDL=nsupercents_by_year_and_country[nsupercents_by_year_and_country$DCOUNTRY=="JPN",]
USA_IDL=nsupercents_by_year_and_country[nsupercents_by_year_and_country$DCOUNTRY=="USA",]

summary(lm(BEL_IDL$n~BEL_IDL$YOD))
#R^2=.004
#p=.33
summary(lm(log(BEL_IDL$n)~BEL_IDL$YOD))
#R^2=.026
#p=.29

summary(lm(ESP_IDL$n~ESP_IDL$YOD))
#R^2=.004
#p=.31
summary(lm(log(ESP_IDL$n)~ESP_IDL$YOD))
#R^2=.003
#p=.32

summary(lm(EW_IDL$n~EW_IDL$YOD))
#R^2=.61
#p=1.6e-9
summary(lm(log(EW_IDL$n)~EW_IDL$YOD))
#R^2=.78
#p=1.8e-14

summary(lm(FRA_IDL$n~FRA_IDL$YOD))
#R^2=.64
#p=1.8e-7
summary(lm(log(FRA_IDL$n)~FRA_IDL$YOD))
#R^2=.78
#p=2.7e-10

summary(lm(JPN_IDL$n~JPN_IDL$YOD))
#R^2=.48
#p=.02
summary(lm(log(JPN_IDL$n)~JPN_IDL$YOD))
#R^2=.73
#p=.02


summary(lm(USA_IDL$n~USA_IDL$YOD))
#R^2=.02
#p=.2
summary(lm(log(USA_IDL$n)~USA_IDL$YOD))
#R^2=-.03
#p=.92


#EW has the best fit and highest sample size
#Although the poor fit in some countries (especially USA!) should give some pause
#EW also spans 1968-2017, which is a good range of years
summary(lm(log(EW_IDL$n)~EW_IDL$YOD))
intercept=-106.44042
coefficient=0.05386
future_years=seq(2021,2100)
predicted_EW_supercents=exp(intercept+future_years*coefficient)

#Assume that EW represents 1% of global supercentenarians
predicted_global_supercents=100*predicted_EW_supercents

cumulative_global_supercents=sum(predicted_global_supercents)
#1,476,796

#Estimates of future benchmarks
prob_115*cumulative_global_supercents
#8267.992
prob_120*cumulative_global_supercents
#40.52564
prob_125*cumulative_global_supercents
#0.2293339

plot(EW_IDL$n~EW_IDL$YOD,ylab="Number of supercentenarians",xlab="Year of Death",xlim=c(1968,2100),pch=19,ylim=c(1,800))
points(predicted_EW_supercents~seq(2021,2100))
legend("bottomright",legend=c("Observed","Projected"),pch=c(19,1))

plot(EW_IDL$n~EW_IDL$YOD,ylab="Number of supercentenarians",xlab="Year of Death",xlim=c(1968,2100),pch=19,ylim=c(1,800),log="y")
points(predicted_EW_supercents~seq(2021,2100))
legend("bottomright",legend=c("Observed","Projected"),pch=c(19,1))
