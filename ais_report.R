#setwd("/Users/mijeongkim/Dropbox/2_research/0_semipara_kernel/R")
library(sn)
library(mixtools)
library(rootSolve)
##########################################################
### functions #############################
reg.m<-function(beta,X){ X%*%beta }
reg.m1<-function(beta,X){ t(X) }

f1f<-function(epsilon,m1,m2,sig1,sig2,p){
  f1<-dnorm(epsilon,m1,sig1)
  f2<-dnorm(epsilon,m2,sig2)
  f1.d<--(epsilon-m1)/sig1^2*f1
  f2.d<--(epsilon-m2)/sig2^2*f2
  
  y<-(p*f1.d+(1-p)*f2.d)/(p*f1+(1-p)*f2)
  return(y)
}

Seff4<-function(theta){
  beta<-theta[1:3]; mu2<-theta[4]
  e<-t(y-reg.m(beta,X))
  n=length(e)
  h=1.06*sd(e)*n^(-1/5);
  fhat=rep(0,n)
  fdhat=rep(0,n)
  for(i in 1:n){
    fhat[i]=sum(dnorm((e[i]-e)/h))/(n*h);
    fdhat[i]=sum(-(e[i]-e)*dnorm((e[i]-e)/h))/(n*h^3);
  }
  mu3<-sum(e^3)/n;
  mu4=sum(e^4)/n;
  t=e^2-mu2-mu3*e/mu2;
  Et2=mu4-mu2^2-mu3^2/mu2;
  f1f=fdhat/fhat;
  
  mat1<-(reg.m1(beta,X)-rowMeans(reg.m1(beta,X)))
  mat2<-rowMeans(reg.m1(beta,X))
  
  Seff_beta<-rbind(-f1f*mat1[1,]+(e/mu2-mu3*t/mu2/Et2)*mat2[1],
                   -f1f*mat1[2,]+(e/mu2-mu3*t/mu2/Et2)*mat2[2],
                   -f1f*mat1[3,]+(e/mu2-mu3*t/mu2/Et2)*mat2[3])
  Seff_mu2<-t/Et2
  rbind(Seff_beta,Seff_mu2)
}

Seff_scaled4<-function(theta){
  rowSums(Seff4(theta))}

estcov4<-function(theta){
  n=length(y)
  p=length(theta)
  f=Seff4(theta)
  B=f%*%t(f)/n
  A=matrix(0,nrow=p,ncol=p)
  delta=1e-6
  for(i in 1:p){
    d=rep(0,p)
    d[i]=theta[i]*delta
    A[,i]=rowMeans(Seff4(theta+d)-Seff4(theta-d))/d[i]/2
  }
  M=diag(solve(A)%*%B%*%t(solve(A)))
  return(M)
}

fe3<-function(m,sig) 3*m*sig^2+m^3
fe4<-function(m,sig) m^4+6*m^2*sig^2+3*sig^4
E.e3<-function(m1,m2,sig1,sig2,p) p*fe3(m1,sig1)+(1-p)*fe3(m2,sig2)
E.e4<-function(m1,m2,sig1,sig2,p) p*fe4(m1,sig1)+(1-p)*fe4(m2,sig2)
E.t2<-function(m1,m2,sig1,sig2,p){
  p*(fe4(m1,sig1)-sig1^2-fe3(m1,sig1)^2/sig1^2)
  +(1-p)*(fe4(m2,sig2)-sig2^2-fe3(m2,sig2)^2/sig2^2)}


Seff<-function(theta,other.param){
  m1<-other.param[1];  m2<-other.param[2]
  sig1<-other.param[3];  sig2<-other.param[4];  p<-other.param[5];
  
  beta<-theta[1:3]; mu2<-theta[4]
  epsilon<-t(y-reg.m(beta,X))
  Ee3<-E.e3(m1,m2,sig1,sig2,p)
  Et2<-E.t2(m1,m2,sig1,sig2,p)
  t<-epsilon^2-mu2-Ee3*epsilon/mu2
  
  mat1<-(reg.m1(beta,X)-rowMeans(reg.m1(beta,X)))
  mat2<-rowMeans(reg.m1(beta,X))
  
  Seff_beta<-rbind(-f1f(epsilon,m1,m2,sig1,sig2,p)*mat1[1,]+(epsilon/mu2-Ee3*t/mu2/Et2)*mat2[1],
                   -f1f(epsilon,m1,m2,sig1,sig2,p)*mat1[2,]+(epsilon/mu2-Ee3*t/mu2/Et2)*mat2[2],
                   -f1f(epsilon,m1,m2,sig1,sig2,p)*mat1[3,]+(epsilon/mu2-Ee3*t/mu2/Et2)*mat2[3])
  Seff_mu2<-t/Et2
  rbind(Seff_beta,Seff_mu2)
}


Seff_scaled<-function(theta,other.param){
  rowSums(Seff(theta,other.param))}

estcov<-function(theta,other.param){
  n=length(y)
  p=length(theta)
  f=Seff(theta,other.param)
  B=f%*%t(f)/n
  A=matrix(0,nrow=p,ncol=p)
  delta=1e-6
  for(i in 1:p){
    d=rep(0,p)
    d[i]=theta[i]*delta
    A[,i]=rowMeans(Seff(theta+d,other.param)-Seff(theta-d,other.param))/d[i]/2
  }
  M=diag(solve(A)%*%B%*%t(solve(A)))
  return(M)
}

################################################################################
### Calculation
### DATA : Australian athelets
data(ais)

### OLSE
Reg1<-lm(Bfat~BMI+LBM,data=ais)
res<-residuals(Reg1)
shapiro.test(res)


#### Semiparametric ####
y<-ais$Bfat
x1<-ais$BMI
x2<-ais$LBM
X<-as.matrix(data.frame(1,x1,x2))

th0<-c(Reg1$coefficients,10)
th4<-multiroot(Seff_scaled4,start=th0)$root  
res4<-y-X%*%th4[1:3]


#### Calculation of SD based on the semiparametric theory
n=length(y)
### Numerical method for calculation of SD
sqrt(estcov4(th4)/n)
### Comparison of mean of residuals
c(mean(res4))
### Comparison of variance of residuals
c(var(res4), th4[4] )

########## Normal mixture ###########
root1<-c(Reg1$coefficients,10)
res.dist<- normalmixEM(res, lambda = .5, mu = c(-3, 3), sigma = 1)
maxiter=20;  
for(iter in 1:maxiter){
  m1<-res.dist$mu[1]; m2<-res.dist$mu[2]
  sig1<-res.dist$sigma[1];  sig2<-res.dist$sigma[2]; p<-res.dist$lambda[1]
  param1<-c(m1,m2,sig1,sig2,p)
  
  root1<-multiroot(Seff_scaled,start=root1,other.param=param1)$root
  res1<-y-X%*%root1[1:3]
  res.dist<- normalmixEM(res1, lambda = p, mu = c(m1, m2), sigma =  c(sig1,sig2))
  #plot(res.dist, density=TRUE)
  #print(iter)
  #print(param1)
  #print(root1)
}
th1<-root1
res1<-y-X%*%th1[1:3]

### Numerical method for calculation of SD
sqrt(estcov(th1,param1)/n)
### Comparison of mean of residuals
c(mean(res1), p*m1+(1-p)*m2)
### Comparison of variance of residuals
c(var(res1), p*(sig1^2+m1^2)+(1-p)*(sig2^2+m2^2), th1[4] )

##########################
### Graph
#pdf("plot2.pdf", width=8.5, height=8)
par(mfrow=c(2,2))

line<-par(lwd=0.1)
hist(res,probability = T,xlim=c(-10,15), xaxt='n',
     xlab="residuals", main=substitute(paste("(A) OLSE
                                             ")));
a<-seq(-10,15,0.01)
axis(side=1,at=seq(-9,15,3))
hist(res4,probability = T,xlim=c(-10,15), xaxt="n",
     xlab="residuals", 
     main=substitute(paste("(B) Semiparametric method with
     kernel density estimation for errors")));
a<-seq(-10,15,0.01)
axis(side=1,at=seq(-9,15,3))
h=1.06*sd(res4)*n^(-1/5);
fhat=rep(0,n)
for(i in 1:n){
  fhat[i]=sum(dnorm((res4[i]-res4)/h))/(n*h);
}
res_dist<-data.frame(res4,fhat)
res_dist1<-res_dist[order(res4),]
lines(res_dist1$res4,res_dist1$fhat,,lwd=1)

### plot for bimodal
hist(res1,breaks=10,probability = T,xlim=c(-10,15), xaxt="n",
     xlab="residuals", 
     main=substitute(paste("(C) Semiparametric method
     with normal mixture errors")));
a<-seq(-10,15,0.01)
ya<-p*dnorm(a,m1,sig1)+(1-p)*dnorm(a,m2,sig2)
axis(side=1,at=seq(-9,15,3))
y1<-p*dnorm(a,m1,sig1);
y2<-(1-p)*dnorm(a,m2,sig2)
lines(a,y1,col="red",lwd=1, lty='twodash')
lines(a,y2,col="blue",lwd=1, lty='twodash')
lines(a,ya, lwd=1)

## QQ plot ###
line<-par(lwd=1)
pn<-((1:n)-0.5)/n
qn<-rep(0,n)
t<-seq(-20,20,0.0001)
pt<-p*pnorm(t,m1,sig1)+(1-p)*pnorm(t,m2,sig2)
for(i in 1:n){
  qn[i]<-t[which.min(abs(pn[i]-pt))]
}
plot(qn,sort(res1),
     ylab='Empirical values', 
     xlab="Theoretical values",
     main=substitute(paste("(D) Q-Q plot of normal mixture errors
                           ")))
lines(c(-10,15),c(-10,15),lty='dashed',col='grey30')

#dev.off()