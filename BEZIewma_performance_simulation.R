#The code below evaluates via Monte Carlo simulation the IC and the OOC performance of a two-sided BEZI EWMA control chart
#The values for the IC design parameters mu0, phi, nu0 as well as the values for lambda and L (values for the chart's design parameters)
#are needed. Then the various RL-based measures are evaluated via simulation.
library(gamlss.dist)
### IC case
phi0<-400 # in-control phi
mu0<-0.01 # in-control mu
nu0<-0.5 # in-control nu
lambda<-0.05 # lambda value
L<-2.473475 # L value, distance of the control limits from center line
#############################
mu0X<-(1-nu0)*mu0
sigma0X2<-(1-nu0)*mu0*(nu0*mu0+(1-mu0)/(1+phi0))
cat("IC mean:",mu0X," IC var:",sigma0X2," IC sd:",sqrt(sigma0X2),"\n") # it prints the values of the BEZI process parameters
########
listRL<-c() # an empty vector
sims<-10000 # number of simulation runs, default 10000
mu0X<-(1-nu0)*mu0
sigma0X2<-(1-nu0)*mu0*(nu0*mu0+(1-mu0)/(1+phi0))
LCL<-mu0X-L*sqrt(sigma0X2*lambda/(2-lambda)) # lower control limit
UCL<-mu0X+L*sqrt(sigma0X2*lambda/(2-lambda)) # upper control limit
####### Simulation Procedure
for(l in 1:sims){
  Z0<-mu0X
  j<-1
  while(TRUE){
    X1<-rBEZI(1, mu =mu0, sigma =phi0, nu =nu0)
    Z1<-lambda*X1+(1-lambda)*Z0
    if(Z1>UCL|Z1<LCL){listRL[l]<-j;break}
    else{
      j<-(j+1)
      Z0<-Z1
    }
  }
}
mean(listRL) # the in-control ARL
sd(listRL) # the in-control SDRL
sd(listRL)/sqrt(sims) # the standard error of the estimated ARL
quantile(listRL,c(0.5,0.95)) # the median run length and the 0.95-percentile point of the run length distribution
######### out-of-control performance -- shifts in mu0 #############
for(delta in c(1.1,1.2,1.3,1.5,1.7,2.0)){
tau<-1 # magnitude of shift in nu0, for tau=1 no shift in nu0
tau2<-0.0 # magnitude of shift in phi0, for tau2=0 no shift in phi0
phi1<-phi0+tau2 # out-of-control phi value
mu1<-delta*mu0 # out-of-control mu value
nu1<-tau*nu0 # out-of-control nu value
#############
mu1X<-(1-nu1)*mu1
sigma1X2<-(1-nu1)*mu1*(nu1*mu1+(1-mu1)/(1+phi1))
cat("IC mean:",mu0X," OOC mean:",mu1X," IC var:",sigma0X2," OOC var:",sigma1X2," IC sd:",sqrt(sigma0X2)," OOC sd:",sqrt(sigma1X2),"\n") # it prints the values of the BEZI process parameters, out-of-control
############### simulation procedure
listRL1<-c() # an empty vector
sims1<-sims # number of simulation runs, same as before
for(l in 1:sims1){
  Z01<-mu0X
  j1<-1
  while(TRUE){
    X11<-rBEZI(1, mu =mu1, sigma =phi0, nu =nu1)
    Z11<-lambda*X11+(1-lambda)*Z01
    if(Z11>UCL|Z11<LCL){listRL1[l]<-j1;break}
    else{
      j1<-(j1+1)
      Z01<-Z11
    }
  }
}
cat(":",delta," :",mean(listRL1)," :",sd(listRL1)," :",sd(listRL1)/sqrt(sims)," :",quantile(listRL1,c(0.5,0.95)),"\n") # performance of the chart
}
########## OOC case -- Shifts in nu0 ######################
for(tau in c(0.9,0.8,0.7,0.5)){
delta<-1.0 # magnitude of shift in mu0, for tau=1 no shift in nu0
tau2<-0.0 # magnitude of shift in phi0, for tau2=0 no shift in phi0
phi1<-phi0+tau2 # out-of-control phi value
mu1<-delta*mu0 # out-of-control mu value
nu1<-tau*nu0 # out-of-control nu value
#############
mu1X<-(1-nu1)*mu1
sigma1X2<-(1-nu1)*mu1*(nu1*mu1+(1-mu1)/(1+phi1))
cat("IC mean:",mu0X," OOC mean:",mu1X," IC var:",sigma0X2," OOC var:",sigma1X2," IC sd:",sqrt(sigma0X2)," OOC sd:",sqrt(sigma1X2),"\n") # it prints the values of the BEZI process parameters, out-of-control
####### simulation procedure
listRL1<-c() # an empty vector
sims1<-sims # number of simulation runs
for(l in 1:sims1){
  Z01<-mu0X
  j1<-1
  while(TRUE){
    X11<-rBEZI(1, mu =mu1, sigma =phi0, nu =nu1)
    Z11<-lambda*X11+(1-lambda)*Z01
    if(Z11>UCL|Z11<LCL){listRL1[l]<-j1;break}
    else{
      j1<-(j1+1)
      Z01<-Z11
    }
  }
}
cat(" :",tau," :",mean(listRL1)," :",sd(listRL1)," :",sd(listRL1)/sqrt(sims)," :",quantile(listRL1,c(0.5,0.95)),"\n") # performance of the chart
}