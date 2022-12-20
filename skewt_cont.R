#setwd("/Users/mijeongkim/Dropbox/2_research/0_semipara_kernel/Matlab/report/skewt_linear_30cont")
library(sn)
library(data.table)
seed=40


set.seed(seed)
#################
simu=1000; n=300;
b0=5; b1=1; b2=1.8;
xi=0; omega=3; alpha=3; nu=10; p=0.7
a=2.5;b=3;
n1=n*p; n2=n*(1-p); 

bv=sqrt(nu)*gamma(0.5*(nu-1))/sqrt(pi)/gamma(0.5*nu);
d=alpha/sqrt(1+alpha^2);
mu=xi+omega*bv*d;
mu2_1=omega^2*(nu/(nu-2)-(bv*d)^2);
mu2_2=a*b^2;
mu2=p*mu2_1+(1-p)*mu2_2;


bhat_set=matrix(NA,simu,3); bcov_set=matrix(NA,simu,3);
bhat_ci<-array(NA,dim=c(simu,2,3)); SEC_para<-matrix(NA,simu,4);
b_olse<-matrix(NA,simu,4)
xy=array(NA,dim=c(n,3,simu)); xy_fails<-NULL
e_set<-matrix(NA,simu,n);
residual_set<-matrix(NA,simu,n);

iter=1; fails=0;
while(iter<=simu){
  bhat=c(NA,NA,NA); bcov=c(NA,NA,NA); check1=NA; m1<-NULL;
  x1<-rgamma(n,5,2); x2<-rgamma(n,4,3);
  e1.1<-rst(n1,xi,omega,alpha,nu); e1<-e1.1-mu;
  e2<-rgamma(n2,a,scale=b)-a*b; e=c(e1,e2);
  y=b0+b1*x1+b2*x2+e; 
  try(m1<- selm(y ~ x1+ x2, family="ST")); 
  try(bhat<-coef(m1,vector=TRUE)[1:3]);
  try(bcov<-diag(vcov(m1))[1:3]); 
  check1<-sum(is.na(bcov))
  if(check1==0){
    bhat_set[iter,]=bhat
    bcov_set[iter,]=bcov
    b_olse[iter,1:3]=coef(reg1<-lm(y~x1+x2)); b_olse[iter,4]=var(residuals(reg1));
    xy[,,iter]<-cbind(x1,x2,y)
    e_set[iter,]<-e
    residual_set[iter,]<-residuals(m1)
    SEC_sum<-summary(extractSECdistr(m1))
    SEC_para[iter,]<-c(-SEC_sum@cp[1],SEC_sum@dp[2:4])
    bhat_ci[iter,,1:3]<-t(confint(m1,parm=1:3))
    iter=iter+1;
    #print(iter)
  }
  else{
    fails=fails+1;
    xy_temp<-cbind(x1,x2,y)
    xy_fails<-array(c(xy_fails,xy_temp),dim=c(n,3,fails))
    #print(fails)
  }
}

bse_set<-sqrt(bcov_set)
bhat_mean<-apply(bhat_set,2,mean)
bhat_median<-apply(bhat_set,2,median)
bhat_sd<-apply(bhat_set,2,sd)
bse_median<-apply(bse_set,2,median)

result<-cbind(bhat_median,bhat_sd,bse_median)

coverage1=sum(between(b0, bhat_ci[,1,1],bhat_ci[,2,1]))/simu*100
coverage2=sum(between(b1, bhat_ci[,1,2],bhat_ci[,2,2]))/simu*100
coverage3=sum(between(b2, bhat_ci[,1,3],bhat_ci[,2,3]))/simu*100

coverage=c(coverage1,coverage2,coverage3)
print(seed)
print(coverage)


save_file_name= paste0("skt1_seed", seed, ".Rdata")
save(xy,xy_fails,result,coverage, b_olse,e_set,residual_set, SEC_para, file = save_file_name)

################################################################################
#### plot #####
#load("skt1_seed40.Rdata")
iter=705

dat1<-xy[,,iter]
x1<-dat1[,1]; x2<-dat1[,2]; y<-dat1[,3]
m1<- selm(y ~ x1+ x2, family="ST")
#plot(m1)
#plot(m1,which=4)

e1<-e_set[iter,]
r1<-residual_set[iter,]
u<-(1:300-0.5)/300
q<-qst(u,SEC_para[iter,1],SEC_para[iter,2],SEC_para[iter,3],SEC_para[iter,4])
v<-c(min(e1)-2,max(c(e1,q))+2)
#plot(sort(e1),q)
#lines(v,v,'l',lty='dashed',col='grey30')

t<-seq(-8,max(e1)+0.5,0.01)
ty1<-dst(t+mu,xi,omega,alpha,nu)
ty2<-dgamma(t+a*b,a,scale=b)
ty=0.7*ty1+0.3*ty2
#plot(t,ty,'l')
#lines(t,ty,'l',col='red',lwd=1)
ty_skt<-dst(t, SEC_para[iter,1],SEC_para[iter,2],SEC_para[iter,3],SEC_para[iter,4])

#######################################################################3
## start drawing plot
#save_plot_name= paste0("skewt_cont30_", iter, ".pdf")
#pdf(save_plot_name, width=8.5, height=8)
pdf("plot1.pdf", width=8.5, height=8)
par(mfrow=c(2,2))

#plot1
line<-par(lwd=0.1)
hist(e1,probability = TRUE,breaks=20,xlim=c(-8,max(e1)+0.5),ylim=c(0,0.20),xaxt='n',
     col=rgb(0.8,0.8,0.8,0.5),xlab="errors",main=substitute(paste("(A) Histogram of errors")))
axis(side=1, at=seq(-8,max(e1)+0.5, 4), labels=seq(-8,max(e1)+0.5, 4))
lines(t,ty,'l',col='red',lwd=1)
lines(t,ty_skt,col='blue',lty='twodash',lwd=1)

#plot2
line<-par(lwd=1)
plot(q,sort(e1),xlim=v,ylim=v,ylab='Empirical values', 
     xlab="Theoretical values",
     main=substitute(paste("(B) Skew-",italic(t)," Q-Q plot of errors")));
lines(v,v,'l',lty='dashed',col='grey30');

#plot3
line<-par(lwd=0.1)
hist(r1,probability = TRUE,breaks=20,xlim=c(-8,max(e1)),ylim=c(0,0.20),xaxt='n',
     col=rgb(0.8,0.8,0.8,0.5),xlab="residuals",main=substitute(paste("(C) Histogram of residuals")))
axis(side=1, at=seq(-8,max(e1), 4), labels=seq(-8,max(e1), 4))
lines(t,ty,'l',col='red',lwd=1)
lines(t,ty_skt,col='blue',lty='twodash',lwd=1)

#plot4
line<-par(lwd=1)
plot(m1,which=4,caption="",
     main=expression(paste("(D) ",chi^2, " P-P plot of (scaled DP residuals)"^2)))

dev.off()






