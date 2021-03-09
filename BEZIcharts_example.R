library(gamlss.dist)
library(gamlss)
library(VGAM)
library(qcc)
### data from Noorossana et al. 2011
# sample sizes ni
n1s<-c(53,45,50,52,45,47,51,52,49,47,49,51,50,52,45,53,48,53,49,53,48,46,49,49,52,49,49,51,46,49,48,50,46,50,53,48,48,50,46,53,47,46,47,46,51,50,47,47,45,50,49,45,48,46,51,45,47,52,50,50,46,49,51,47,50,51,48,49,45,53,49,51,49,49,50,46,50,48,49,51,48,49,52,45,47,46,48,51,45,50,50,48,52,47,49,50,53,53,52,46)
### number of daily infections Xi
x0<-c(2,      0,      0,      0,      1,      0,      2,      0,      4,      0,      0,      0,      0,      0,      0,      0,
      2,      0,      3,      0,      1,      0,      0,      0,      0,      0,      1,      1,      0,      0,      0,      3,
      0,      0,      4,      0,      2,      0,      0,      2,      0,      4,      1,      0,      0,      3,      0,      1,
      0,      4,      0,      0,      3,      0,      0,      0,      0,      0,      1,      0,      0,      0,      1,      3,
      0,      0,      0,      5,      0,      0,      1,      0,      0,      1,      0,      0,      2,      2,      1,      1,
      3,      0,      0,      1,      2,      0,      0,      0,      0,      0,      3,      0,      0,      0,      0,      0,
      0,      0,      1,      0)
p0j<-x0/n1s # daily proportions, Wi values
######## a p-chart with 3sigma limits, using the package qcc
pchartpic<-qcc(x0,type="p",sizes=n1s,plot=FALSE)
plot(pchartpic,add.stats=TRUE,title="p-chart",xlab="sample",ylab="proportion")
### fit bezi distribution on the available data, with the function gamlssML (maximum likelihood estimation)
fitbezi<-gamlssML(p0j,family=BEZI)
summary(fitbezi)
####
logitmu<-coef(fitbezi,what="mu")
muest<-exp(logitmu)/(1+exp(logitmu))
logsigma<-coef(fitbezi,what="sigma")
sigmaest<-exp(logsigma)
logitnu<-coef(fitbezi,what="nu")
nuest<-exp(logitnu)/(1+exp(logitnu))
cat("mu:",muest," phi:",sigmaest," nu:",nuest,"\n") # estimated parameters mu, phi and nu
##################### BEZI-Shewhart, control chart construction
plot(1:length(p0j),p0j,ylim=c(0,0.16),type="n",main="BEZI-Shewhart chart",xlab="sample",ylab="proportion")
points(1:length(p0j),p0j,pch=18)
lines(1:length(p0j),p0j,lty=1)
abline(h=0.01466427,lty=2) # center line at mu0
abline(h=0.1157912,lty=1) # upper control limit, ARL0=370.
text(90,0.125,"UCL=0.1157912")
####### BEZI-EWMA, control chart construction
par(mfrow=c(2,2)) # four graphs in 2 columns and 2 rows
############# 
lam1<-0.05 # lambda=0.05
Lewma<-2.466107 # L value
mu=0.0431302 # in-control mu
phi=77.6967 # in-control phi
nu=0.66 # in-control nu
mu0X<-(1-nu)*mu # in-control process mean level
sigma0X2<-(1-nu)*mu*(nu*mu+(1-mu)/(1+phi)) # in-control variance
CLewma<-mu0X # center line
UCLewma<-mu0X+Lewma*sqrt(sigma0X2*lam1/(2-lam1)) # upper control limit, EWMA chart
LCLewma<-mu0X-Lewma*sqrt(sigma0X2*lam1/(2-lam1)) # lower control limit, EWMA chart
z0<-mu0X # starting value
listz<-c() # empty vector
for(i in 1:length(p0j)){
z1<-lam1*p0j[i]+(1-lam1)*z0
listz[i]<-z1
z0<-z1
}
## the for calculates the values of the EWMA statisti
plot(1:length(p0j),listz,type="n",
     main=expression("BEZI-EWMA chart,"~lambda==0.05),
     xlab="sample",
     ylab=expression("EWMA Z"[i]),
     ylim=c(0,0.035))
points(1:length(p0j),listz,pch=18)
lines(1:length(p0j),listz,lty=1)
abline(h=UCLewma,lty=1)
abline(h=max(0,LCLewma),lty=1)
abline(h=CLewma,lty=2)
text(15,0.0275,"UCL=0.02430269")
text(15,0.0025,"LCL=0.005025847")
which(listz>UCLewma) # which points give an OOC signal
print(c(LCLewma,CLewma,UCLewma)) # prints the UCL, LCL and CL values
####################################### 
## below we construct the BEZI-EWMA charts for lambda = 0.10, 0.20, 0.30
## the code is the same as for the case of lambda = 0.05
## changes are needed in lambda, L values as well as for the control limits.
lam1<-0.10
Lewma<-2.82306
mu=0.0431302
phi=77.6967
nu=0.66
mu0X<-(1-nu)*mu
sigma0X2<-(1-nu)*mu*(nu*mu+(1-mu)/(1+phi))
CLewma<-mu0X
UCLewma<-mu0X+Lewma*sqrt(sigma0X2*lam1/(2-lam1))
LCLewma<-mu0X-Lewma*sqrt(sigma0X2*lam1/(2-lam1))
z0<-mu0X
listz<-c()
for(i in 1:length(p0j)){
z1<-lam1*p0j[i]+(1-lam1)*z0
listz[i]<-z1
z0<-z1
}
plot(1:length(p0j),listz,type="n",
     main=expression("BEZI-EWMA chart,"~lambda==0.10),
     xlab="sample",
     ylab=expression("EWMA Z"[i]),
     ylim=c(0,0.036))
points(1:length(p0j),listz,pch=18)
lines(1:length(p0j),listz,lty=1)
abline(h=UCLewma,lty=1)
abline(h=max(0,LCLewma),lty=1)
abline(h=CLewma,lty=2)
text(9,0.033,"UCL=0.03047")
which(listz>UCLewma)
print(c(LCLewma,CLewma,UCLewma))
##########################
lam1<-0.20
Lewma<-3.275712
mu=0.0431302
phi=77.6967
nu=0.66
mu0X<-(1-nu)*mu
sigma0X2<-(1-nu)*mu*(nu*mu+(1-mu)/(1+phi))
CLewma<-mu0X
UCLewma<-mu0X+Lewma*sqrt(sigma0X2*lam1/(2-lam1))
LCLewma<-mu0X-Lewma*sqrt(sigma0X2*lam1/(2-lam1))
z0<-mu0X
listz<-c()
for(i in 1:length(p0j)){
z1<-lam1*p0j[i]+(1-lam1)*z0
listz[i]<-z1
z0<-z1
}
plot(1:length(p0j),listz,type="n",
     main=expression("BEZI-EWMA chart,"~lambda==0.20),
     xlab="sample",
     ylab=expression("EWMA Z"[i]),
     ylim=c(0,0.05))
points(1:length(p0j),listz,pch=18)
lines(1:length(p0j),listz,lty=1)
abline(h=UCLewma,lty=1)
abline(h=max(0,LCLewma),lty=1)
abline(h=CLewma,lty=2)
text(12,0.046,"UCL=0.04132")
which(listz>UCLewma)
print(c(LCLewma,CLewma,UCLewma))
##################
lam1<-0.30
Lewma<-3.53015
mu=0.0431302
phi=77.6967
nu=0.66
mu0X<-(1-nu)*mu
sigma0X2<-(1-nu)*mu*(nu*mu+(1-mu)/(1+phi))
CLewma<-mu0X
UCLewma<-mu0X+Lewma*sqrt(sigma0X2*lam1/(2-lam1))
LCLewma<-mu0X-Lewma*sqrt(sigma0X2*lam1/(2-lam1))
z0<-mu0X
listz<-c()
for(i in 1:length(p0j)){
z1<-lam1*p0j[i]+(1-lam1)*z0
listz[i]<-z1
z0<-z1
}
plot(1:length(p0j),listz,type="n",
     main=expression("BEZI-EWMA chart,"~lambda==0.30),
     xlab="sample",
     ylab=expression("EWMA Z"[i]),
     ylim=c(0,0.07))
points(1:length(p0j),listz,pch=18)
lines(1:length(p0j),listz,lty=1)
abline(h=UCLewma,lty=1)
abline(h=max(0,LCLewma),lty=1)
abline(h=CLewma,lty=2)
text(9,0.0558,"UCL=0.05086")
which(listz>UCLewma)
print(c(LCLewma,CLewma,UCLewma))
##########################
par(mfrow=c(1,1))