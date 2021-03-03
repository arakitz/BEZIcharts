#The code below evaluates the performance of the two-sided BEZI EWMA control chart
#The values of the process parameters mu0, nu0 and phi must be given as well as the
#values of the chart's design parameters lambda and L
#The Markov chain method is used for evaluating the performance of the chart
library(gamlss.dist)
library(nleqslv)
#library(FAdist)
############### parameters
m<-200 # states of the Markov chain
M<-2*m+1
phi0<-400 # in-control value for phi
mu0<-0.01 # in-control value for parameter mu
nu0<-0.5 # in-control value for the probability of zero occurrence
plotBEZI(mu=mu0,sigma=phi0,nu=nu0,from=0,to=0.999,n=101) # a plot of the IC BEZI distribution
######################
lambda<-0.05 # parameter lambda for the EWMA chart
L<-2.473475 # parameter L for the EWMA chart, distance of the control limits from the center line
##### calculcation of control limits
mu0X<-(1-nu0)*mu0
sigma0X2<-(1-nu0)*mu0*(nu0*mu0+(1-mu0)/(1+phi0))
LCL<-mu0X-L*sqrt(sigma0X2*lambda/(2-lambda))
UCL<-mu0X+L*sqrt(sigma0X2*lambda/(2-lambda))
########################### Markov chain Method #########################
delta<-(UCL-LCL)/(M)
Q0<-matrix(0,ncol=M,nrow=M)
  for(i in 1:M){
    for(j in 1:M){
      lowerl<-LCL+delta*(j-1-(1-lambda)*(i-0.5))/lambda
      upperl<-LCL+delta*(j-(1-lambda)*(i-0.5))/lambda
      Q0[i,j]<-pBEZI(upperl,mu = mu0, sigma = phi0, nu = nu0)-pBEZI(lowerl,mu = mu0, sigma = phi0, nu = nu0)
    }
  }
  avec<-rep(0,M);avec[m+1]<-1
  l1<-rep(1,M)
  ID<-diag(M)
  ARLin<-avec%*%solve(ID-Q0)%*%l1
#### out-of-control performance, changes in nu0 #############  
for(dl in c(1.0,0.9,0.8,0.7,0.5)){
  mu1<-1.0*mu0 # no change in mu0
  phi1<-1.0*phi0 # no change in phi0
  nu1<-dl*nu0 # the out-of-control value for nu
####
  mu1X<-(1-nu1)*mu1
  sigma1X2<-(1-nu1)*mu1*(nu1*mu1+(1-mu1)/(1+phi1))
####
  Q1u<-matrix(0,ncol=M,nrow=M)
  for(i in 1:M){
    for(j in 1:M){
      lowerlu<-LCL+delta*(j-1-(1-lambda)*(i-0.5))/lambda
      upperlu<-LCL+delta*(j-(1-lambda)*(i-0.5))/lambda
      Q1u[i,j]<-pBEZI(upperlu,mu = mu1, sigma = phi1, nu = nu1)-pBEZI(lowerlu,mu = mu1, sigma = phi1, nu = nu1)
    }
  }
  ARLout<-avec%*%solve(ID-Q1u)%*%l1
  Mout<-solve(ID-Q1u)
  E2RLout<-ARLout+(2*avec%*%Mout%*%Mout%*%Q1u%*%l1)
  SDRLout<-sqrt(E2RLout-ARLout^2)
  m1out<-ARLout
  m2out<-E2RLout
  ERL3out<-6*avec%*%Mout%*%Mout%*%Mout%*%Q1u%*%Q1u%*%l1
  m3out<-ERL3out+3*E2RLout-2*m1out
  MM1<-m1out
  MM2<-m2out-m1out^2
  MM3<-m3out-3*m2out*m1out+2*(m1out^3)
  A1<-4*(MM2^3)/(MM3^2)
  B1<-MM3/(2*MM2)
  threspar<-MM1-(2*MM2^2)/(MM3)
  MRL<-qgamma(0.5,shape=A1,scale=B1)+threspar # calculation of the median run length, MRL
  RL95<-qgamma(0.95,shape=A1,scale=B1)+threspar # calculation of the 0.95-percentile point of the run-length distribution
  ##### output for the results
cat("L:",L," shift:",dl," ARL0:",ARLin," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",ceiling(MRL)," RL095:",ceiling(RL95)," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu1:",mu1X," varX1:",sigma1X2,"\n")
}
#### out-of-control performance, changes in mu0 #############
for(ds in c(1.0,1.1,1.2,1.3,1.5,1.7,2.0)){
  mu1<-ds*mu0 # out-of-control value for mu
  phi1<-1.0*phi0 # no change in phi
  nu1<-1.0*nu0 # no change in nu
####
  mu1X<-(1-nu1)*mu1
  sigma1X2<-(1-nu1)*mu1*(nu1*mu1+(1-mu1)/(1+phi1))
####
  Q1u<-matrix(0,ncol=M,nrow=M)
  for(i in 1:M){
    for(j in 1:M){
      lowerlu<-LCL+delta*(j-1-(1-lambda)*(i-0.5))/lambda
      upperlu<-LCL+delta*(j-(1-lambda)*(i-0.5))/lambda
      Q1u[i,j]<-pBEZI(upperlu,mu = mu1, sigma = phi1, nu = nu1)-pBEZI(lowerlu,mu = mu1, sigma = phi1, nu = nu1)
    }
  }
  ARLout<-avec%*%solve(ID-Q1u)%*%l1
  Mout<-solve(ID-Q1u)
  E2RLout<-ARLout+(2*avec%*%Mout%*%Mout%*%Q1u%*%l1)
  SDRLout<-sqrt(E2RLout-ARLout^2)
  m1out<-ARLout # E(RL)
  m2out<-E2RLout # E(RL^2)
  ERL3out<-6*avec%*%Mout%*%Mout%*%Mout%*%Q1u%*%Q1u%*%l1
  m3out<-ERL3out+3*m2out-2*m1out # E(RL^3)
  MM1<-m1out # mu1
  MM2<-m2out-(m1out^2) # mu2
  MM3<-m3out-3*m2out*m1out+2*(m1out^3) # mu3
  A1<-4*(MM2^3)/(MM3^2)
  B1<-MM3/(2*MM2)
  threspar<-MM1-(2*MM2^2)/(MM3)
  MRL<-qgamma(0.5,shape=A1,scale=B1)+threspar # calculation of the median run length, MRL
   RL95<-qgamma(0.95,shape=A1,scale=B1)+threspar # calculation of the 0.95-percentile point of the run-length distribution
  ##### output of the results #############
  cat("L:",L," shift:",ds," ARL0:",ARLin," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",ceiling(MRL)," RL095:",ceiling(RL95)," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu1:",mu1X," varX1:",sigma1X2,"\n")
}

