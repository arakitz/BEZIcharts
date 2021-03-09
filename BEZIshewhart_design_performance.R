library(gamlss.dist)
library(nleqslv)
############### parameters of the BEZI distribution
phi0<-400 # IC phi value
mu0<-0.1 # IC mu value
nu0<-0.8 # IC nu value
plotBEZI(mu=mu0,sigma=phi0,nu=nu0,from=0,to=0.999,n=101) # plot of the pdf of the BEZI distribution
######################
mu0X<-(1-nu0)*mu0 # IC mean of the BEZI process
sigma0X2<-(1-nu0)*mu0*(nu0*mu0+(1-mu0)/(1+phi0)) # IC variance of the BEZI process
######## Determination of the control limits of the BEZI Shewhart chart
LCL0<-qBEZI(0.00135,mu=mu0,sigma=phi0,nu=nu0) # lower control limit, only in the case of alpha/2>nu0
LCL<-ifelse(LCL0>0,LCL0,LCL<-(LCL0-1))
# below we determine the upper control limit of the BEZI-Shewhart chart
UCL<-ifelse(LCL0>0,UCL<-qBEZI(1-0.00135,mu=mu0,sigma=phi0,nu=nu0),UCL<-qBEZI(1-0.0027,mu=mu0,sigma=phi0,nu=nu0))
alpha<-1-(pBEZI(UCL,mu=mu0,sigma=phi0,nu=nu0)-pBEZI(LCL,mu=mu0,sigma=phi0,nu=nu0)) # the FAR
alpha # prints the FAR
##################### OOC shifts in nu, performance of the BEZI-Shewhart chart
for(dl in c(1.0,0.9,0.8,0.7,0.5)){
mu1<-1.0*mu0 # OOC mu value, same as the IC value
phi1<-phi0 # OOC phi value, same as the IC value
nu1<-dl*nu0 # OOC nu value,
####
mu1X<-(1-nu1)*mu1 # OOC mean of the BEZI process
sigma1X2<-(1-nu1)*mu1*(nu1*mu1+(1-mu1)/(1+phi1)) # OOC variance of the BEZI process
####
beta1<-1-(pBEZI(UCL,mu=mu1,sigma=phi1,nu=nu1)-pBEZI(LCL,mu=mu1,sigma=phi1,nu=nu1)) # OOC probability for point beyond control limit(s)
ARLout<-1/beta1 # OOC ARL for the BEZI-Shewhart chart
SDRLout<-sqrt(1-beta1)/beta1 # OOC SDRL for the BEZI-Shewhart chart
MRL<-ceiling(log(1-0.5)/log(1-beta1)) # OOC MRL for the BEZI-Shewhart chart
RL95<-ceiling(log(1-0.95)/log(1-beta1)) # OOC 0.95-percentile point for the BEZI-Shewhart chart
  cat(" shift:",dl," ARL0:",1/alpha," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",MRL," RL95:",RL95," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu1:",mu1X," varX1:",sigma1X2,"\n")
}
######### OOC shifts in mu, performance of the BEZI-Shewhart chart
for(ds in c(1.0,1.1,1.2,1.3,1.5,1.7,2.0)){
  mu1<-ds*mu0 # OOC mu value
  phi1<-phi0 # OOC phi value, same as the IC value
  nu1<-1.0*nu0 # OOC nu value, same as the IC value
####
  mu1X<-(1-nu1)*mu1 # OOC mean of the BEZI process
  sigma1X2<-(1-nu1)*mu1*(nu1*mu1+(1-mu1)/(1+phi1)) # OOC variance of the BEZI process
####
  beta1<-1-(pBEZI(UCL,mu=mu1,sigma=phi1,nu=nu1)-pBEZI(LCL,mu=mu1,sigma=phi1,nu=nu1)) # OOC probability for point beyond control limit(s)
  ARLout<-1/beta1 # OOC ARL for the BEZI-Shewhart chart
SDRLout<-sqrt(1-beta1)/beta1 # OOC SDRL for the BEZI-Shewhart chart
MRL<-ceiling(log(1-0.5)/log(1-beta1)) # OOC MRL for the BEZI-Shewhart chart
RL95<-ceiling(log(1-0.95)/log(1-beta1)) # OOC 0.95-percentile point for the BEZI-Shewhart chart
  #  ARLout
  cat(" shift:",ds," ARL0:",1/alpha," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",MRL," RL95:",RL95," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu1:",mu1X," varX1:",sigma1X2,"\n")
}
