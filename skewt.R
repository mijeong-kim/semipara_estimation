#setwd("/Users/mijeongkim/Dropbox/2_research/0_semipara_kernel/Matlab/report/skewt_linear")
library(sn)
library(data.table)
seed=75


set.seed(seed)
#################
simu=1000; n=300;
b0=5; b1=1; b2=1.8;
xi=0; omega=3; alpha=3; nu=10;
b_olse<-matrix(NA,simu,4);
bhat_set=matrix(NA,simu,3); bcov_set=matrix(NA,simu,3); 
bhat_ci<-array(NA,dim=c(simu,2,3)); SEC_para<-matrix(NA,simu,4);
xy=array(NA,dim=c(n,3,simu));  xy_fails<-NULL;
bv=sqrt(nu)*gamma(0.5*(nu-1))/sqrt(pi)/gamma(0.5*nu);
d=alpha/sqrt(1+alpha^2);
mu=xi+omega*bv*d; 

iter=1; n_MPLE=0; fails=0;
while(iter<=simu){
  bhat=c(NA,NA,NA); bcov=c(NA,NA,NA); check1=NA; m1<-NULL;
  x1<-rgamma(n,5,2); x2<-rgamma(n,4,3);
  e1<-rst(n,xi,omega,alpha,nu); e<-e1-mu;
  y=b0+b1*x1+b2*x2+e; 
  try(m1<- selm(y ~ x1+ x2, family="ST")); 
  try(bhat<-coef(m1,vector=TRUE)[1:3]);
  try(bcov<-diag(vcov(m1))[1:3]); 
  check1<-sum(is.na(bcov))
  if(check1==0){
    #print(iter)
    bhat_set[iter,]=bhat
    bcov_set[iter,]=bcov
    b_olse[iter,1:3]=coef(reg1<-lm(y~x1+x2)); b_olse[iter,4]=var(residuals(reg1));
    xy[,,iter]<-cbind(x1,x2,y)
    SEC_sum<-summary(extractSECdistr(m1))
    SEC_para[iter,]<-c(-SEC_sum@cp[1],SEC_sum@dp[2:4])
    bhat_ci[iter,,1:3]<-t(confint(m1,parm=1:3))
    iter=iter+1;
    }else{
        fails=fails+1;
        #print(fails)
        xy_temp<-cbind(x1,x2,y)
        xy_fails<-array(c(xy_fails,xy_temp),dim=c(n,3,fails)) 
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
print(coverage)


######
save_file_name= paste0("seed", seed, ".Rdata")
save(xy,result,coverage,b_olse, SEC_para, fails, xy_fails,file = save_file_name)

