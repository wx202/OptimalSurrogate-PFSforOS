rm(list=ls())
timestart<-Sys.time()
# setwd("/Users/celehs_hms_desktop/Dropbox (Harvard University)/Xuan/Surrogate_BothCensored_ProgressionFreeSurvival2023/simu")
setwd("/Users/xuanwang/Dropbox (Harvard University)/Xuan/Surrogate_BothCensored_PFS/simu")
# setwd("/home/xw127/")
source("funs02052024-06292025-github.R")
library(survival)
library(SurrogateOutcome)

# data=data.frame(xob,deltaob,aob,sob)
# save(data,file='data.rda')
load("~/Harvard University Dropbox/Wang Xuan/Xuan/Surrogate_BothCensored_PFS/simu/data.rda")

t.0=1; t=5; stept=0.1; tt=seq(t.0,t,stept)
nn=200; re=5

out=est(t,t.0,tt,nn=200,re=100,data)


timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime) 
