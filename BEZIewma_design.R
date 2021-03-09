library(gamlss.dist)
library(nleqslv)
############### parameters of the BEZI process
phi0<-77.6967 # in-control value of parameter phi
mu0<-0.0431302 # in-control value of parameter mu
nu0<-0.66 # in-control value for the zero-inflated parameter
plotBEZI(mu=mu0,sigma=phi0,nu=nu0,from=0,to=0.999,n=101) # pdf of the BEZI distribution
###################### Design state, Markov chain method
m<-200 # states for the Markov chain
M<-2*m+1 # matrix dimension, number of subintervals in the interval [LCL,UCL]
# for each value lambda, we use the nleqslv function in R to find the unique value for L
# we define a function, named myfun2, for this purpose
# the function is for the case ARL0=370.4. It needs a slight modification to work for any other ARL0 value
for(lambda in seq(0.30,0.30,by=0.1)){
  myfun2<-function(y){
    #### below we determine the entries of the Transition Probability Matrix Q  
    mu0X<-(1-nu0)*mu0
    sigma0X2<-(1-nu0)*mu0*(nu0*mu0+(1-mu0)/(1+phi0))
    LCLy<-mu0X-y*sqrt(sigma0X2*lambda/(2-lambda))
    UCLy<-mu0X+y*sqrt(sigma0X2*lambda/(2-lambda))
    deltay<-(UCLy-LCLy)/(M)
    Q0y<-matrix(0,ncol=M,nrow=M)
    for(i in 1:M){
      for(j in 1:M){
        lowerly<-LCLy+deltay*(j-1-(1-lambda)*(i-0.5))/lambda
        upperly<-LCLy+deltay*(j-(1-lambda)*(i-0.5))/lambda
        Q0y[i,j]<-pBEZI(upperly,mu = mu0, sigma = phi0, nu = nu0)-pBEZI(lowerly,mu = mu0, sigma = phi0, nu = nu0)
      }
    }
    avecy<-rep(0,M);avecy[m+1]<-1 # initial probabilities vector
    l1y<-rep(1,M) # vector of ones
    IDy<-diag(M) # Identity Matrix
    ARLyin<-avecy%*%solve(IDy-Q0y)%*%l1y # ARL
    ARLyin-370.4 # we define ARL-ARL0 and with nleqslv we solve the equation ARL-ARL0=0
  }
  obj2<-nleqslv(c(2.81),myfun2,method="Newton",global="qline",control=list(btol=0.00001,maxit=100)) # solution of the equation ARL-ARL0=0
  L<-obj2$x # the L value
######## IN Control
# verification, for the L value, the ARL0 should be 370.4
  mu0X<-(1-nu0)*mu0 # IC mean of the BEZI distribution
  sigma0X2<-(1-nu0)*mu0*(nu0*mu0+(1-mu0)/(1+phi0)) # IC variance of the BEZI distribution
  LCL<-mu0X-L*sqrt(sigma0X2*lambda/(2-lambda))
  UCL<-mu0X+L*sqrt(sigma0X2*lambda/(2-lambda))
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
  print(c(lambda,L,as.vector(ARLin))) # here the output is lambda, L and the ARL0 value
# if everything is OK, the ARL0 should be 370.4  
}
#}
######### OOC
#### we evaluate the OOC performance of the BEZI-EWMA chart for various shifts in process parameters
LCL<-mu0X-L*sqrt(sigma0X2*lambda/(2-lambda))
UCL<-mu0X+L*sqrt(sigma0X2*lambda/(2-lambda))
delta<-(UCL-LCL)/(M)
### Below is the case of shifts only in nu0
for(dl in c(1.0,0.9,0.8,0.7,0.5)){
  mu1<-1.0*mu0 # mu1 value, same as the IC value mu0
  phi1<-1.0*phi0 # phi1 value, same as the IC value phi0
  nu1<-dl*nu0 # nu1 OOC value
  ####
  mu1X<-(1-nu1)*mu1 # OOC mean of the BEZI distribution
  sigma1X2<-(1-nu1)*mu1*(nu1*mu1+(1-mu1)/(1+phi1)) # OOC variance of the BEZI distribution
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
  MRL<-qgamma(0.5,shape=A1,scale=B1)+threspar # the MRL value
  RL95<-qgamma(0.95,shape=A1,scale=B1)+threspar # the 0.95-percentile value of the run-length distribution
#####
  cat("L:",L," shift:",dl," ARL0:",ARLin," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",ceiling(MRL)," RL095:",ceiling(RL95)," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu1:",mu1X," varX1:",sigma1X2,"\n")
}
### shift in mu
for(ds in c(1.0,1.1,1.2,1.3,1.5,1.7,2.0)){
  mu1<-ds*mu0 # mu1 value, OOC
  phi1<-1.0*phi0 # phi1 value, it is the same as the IC value phi0
  nu1<-1.0*nu0 # nu1 values, it is the same as the IC value nu0
  ####
  mu1X<-(1-nu1)*mu1 # OOC mean of the BEZI distribution
  sigma1X2<-(1-nu1)*mu1*(nu1*mu1+(1-mu1)/(1+phi1)) # OOC variance of the BEZI distribution
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
  MRL<-qgamma(0.5,shape=A1,scale=B1)+threspar # the MRL value
   RL95<-qgamma(0.95,shape=A1,scale=B1)+threspar # the 0.95-percentile value of the run-length distribution
  #####
  cat("L:",L," shift:",ds," ARL0:",ARLin," ARL1:",ARLout," SDRL1:",SDRLout," MRL:",ceiling(MRL)," RL095:",ceiling(RL95)," LCL:",round(LCL,digits=5)," UCL:",round(UCL,digits=5)
      ," mu1:",mu1X," varX1:",sigma1X2,"\n")
}