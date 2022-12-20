function table1_n500()
for p=29
randn('state',p); rand('state',p);
simus=1000;n=500; s=0.9;
beta0=[12 -0.5]';mu2=s^2*pi^2/3;
theta0=[beta0; mu2];
th1=zeros(2,simus); 
th2=zeros(3,simus); th3=zeros(3,simus); th4=zeros(3,simus);
cov2=zeros(3,simus); cov3=zeros(3,simus);  cov4=zeros(3,simus);
samplevar=zeros(4,simus);
e1_set=zeros(n,simus);  e2_set=zeros(n,simus);  e3_set=zeros(n,simus);  e4_set=zeros(n,simus);
H=zeros(1,simus); pValue=zeros(1,simus); W=zeros(1,simus);

%options2=optimset('MaxFunEvals',1e4,'display','off');
options3=optimset('algorithm','levenberg-marquardt','MaxFunEvals',1e4,'MaxIter',1e4,'display','off');
%options3=optimset('MaxFunEvals',1e4,'MaxIter',1e4,'display','off');
p
iter=1;
fails=0;
while (iter<=simus)
  [x,y]=gendata(s,beta0,n);
 eflag=1;
 %OLS
 if (eflag>0)
   [th1(:,iter),f,eflag]=fminsearch(@e2sum,beta0,[],x,y);
    e1=y-m(x,th1(:,iter));
    if (eflag<=0) 
     [iter 1 eflag th1(:,iter)']%'
     fails=fails+1;
   end
 end
%STARTING Semiparametric estimator
 if (eflag>0)
     [th2(:,iter),f,eflag]=fsolve(@seffscaled,theta0,options3,x,y);
     cov2(:,iter)=estcov_d(th2(:,iter),x,y);
     e2=y-m(x,th2(:,iter));
     if (eflag<=0)
       eflag=0;
       fails=fails+1;
     end
 end
  if (eflag>0)
     [th3(:,iter),f,eflag]=fsolve(@seffscaled3,theta0,options3,x,y);
     cov3(:,iter)=estcov3(th3(:,iter),x,y);
     e3=y-m(x,th3(:,iter));
     if ((eflag<=0))
         %|(cov3(3,iter)>200))
         %[iter 6 eflag cov3(:,iter)']%'
       eflag=0;
       fails=fails+1;
     end
  end
    if (eflag>0)
     [th4(:,iter),f,eflag]=fsolve(@seffscaled4,theta0,options3,x,y);
     cov4(:,iter)=estcov4_d(th4(:,iter),x,y);
     e4=y-m(x,th4(:,iter));
     if (eflag<=0)
       eflag=0;
       fails=fails+1;
     end
 end
 if (eflag>0)
    e1_set(:,iter)=e1;      samplevar(1,iter)=var(e1);
    e2_set(:,iter)=e1;      samplevar(2,iter)=var(e2);
    e3_set(:,iter)=e1;      samplevar(3,iter)=var(e3);  [H(iter),pValue(iter),W(iter)]=swtest(e3,0.05);
    e4_set(:,iter)=e1;      samplevar(4,iter)=var(e4);      
  iter=iter+1
 end
end
fails

coverage2=cov_cal(th2,cov2,theta0)/simus*100
coverage3=cov_cal(th3,cov3,theta0)/simus*100
coverage4=cov_cal(th4,cov4,theta0)/simus*100


[median(th2,2) std(th2')'   median(sqrt(cov2),2)]
[median(th3,2) std(th3')'   median(sqrt(cov3),2)]
[median(th4,2) std(th4')'   median(sqrt(cov4),2)]
save table1_n500.mat th1 th2 th3  th4 cov2 cov3  cov4 samplevar p s n coverage2 coverage3 coverage4 e1_set e2_set e3_set e4_set e4_set H pValue W
end

function [x,y]=gendata(s,beta,n)
%x=unifrnd(0,10,1,n);
%x=chi2rnd(5,1,n);
pd1= makedist('Gamma','a',2.5,'b',1.5);
x=random(pd1,1,n);
pd = makedist('Logistic','mu',0,'sigma',s);
e=random(pd,1,n);     
y=m(x,beta)+e;

function f=e2sum(th,x,y)
f=sum((y-m(x,th)).^2);

function f=f1f(e,s)
f=-1/s+2/s*1./(1+exp(-e./s));


function y=m(x,beta)
y=beta(1)*exp(beta(2)*x);

function y=dm(x,theta)
y=[exp(theta(2)*x) 
    m(x,theta).*x];

function f=seff(th,x,y)
e=y-m(x,th); 
%pd=fitdist(e','Logistic');
%mu=pd.mu; s=pd.sigma;
options3=optimset('MaxFunEvals',1e4,'MaxIter',1e4,'display','off');
%s=fminsearch(@loglogistic,0.9,options3,e);
mu2=th(3); mu3=0;
mu4=6/5*mu2^2;
%mu4=mean(e.^4);
s=sqrt(3*mu2)/pi;
t=e.^2-mu2-mu3.*e./mu2;
Et2=mu4-mu2^2-mu3^2/mu2;
f=[-f1f(e,s).*(dm(x,th)-mean(dm(x,th),2))+(e./mu2-mu3*t./(mu2*Et2)).*mean(dm(x,th),2);
    t./Et2];

function f=seff3(th,x,y)
e=y-m(x,th);  
mu2=th(3);
t2=2*mu2.^2;
f=[e./mu2.*dm(x,th);
    (e.^2-mu2)./t2];

function f=loglogistic(s,e)
f=-sum(-e./s-log(s*(1+exp(-e./s)).^2));
 

function f=seff4(th,x,y)
e=y-m(x,th); 
n=size(e,2);
%h=0.9*min(std(e),iqr(e)/1.34)*n^(-1/5);
h=1.06*std(e)*n^(-1/5);
mu2=th(3);
fhat=zeros(size(e));
fdhat=zeros(size(e));
for i=1:n
    fhat(i)=sum(normpdf((e(i)-e)./h))./(n*h);
    fdhat(i)=sum(-(e(i)-e).*normpdf((e(i)-e)./h))./(n*h^3);
end
%mu3=sum(e.^3+3*e*h^2)/n;
mu3=sum(e.^3)/n;
%mu4=sum(e.^4+6*e.^2*h^2+3*h^4)/n;
mu4=sum(e.^4)/n;
t=e.^2-mu2-mu3.*e./mu2;
Et2=mu4-mu2.^2-mu3.^2./mu2;
f1f=fdhat./fhat;
f=[-f1f.*(dm(x,th)-mean(dm(x,th),2))+(e./mu2-mu3*t./mu2./Et2).*mean(dm(x,th),2);
    t./Et2];


function f=seffscaled(th,x,y)
f=mean(seff(th,x,y),2);


function f=seffscaled3(th,x,y)
f=mean(seff3(th,x,y),2);

function f=seffscaled4(th,x,y)
f=mean(seff4(th,x,y),2);

function cov=estcov_d(th,x,y)
n=length(y);
p=length(th);
f=seff(th,x,y);
B=f*f'/n;%'  
delta=1e-6;
for i=1:p
 d=zeros(p,1);					
 d(i)=th(i)*delta;
 A(:,i)=mean(seff(th+d,x,y)-seff(th-d,x,y),2)/d(i)/2;
end
cov=diag((A\B)/(A')/n);%'

function cov=estcov4(th,x,y)
n=length(y);
p=length(th);
f=seff4(th,x,y);
B=f*f'/n;%' 
cov=diag(inv(B)/n);

function cov=estcov3(th,x,y)
n=length(y);
p=length(th);
f=seff3(th,x,y);
B=f*f'/n;%'  
delta=1e-6;
for i=1:p
 d=zeros(p,1);					
 d(i)=th(i)*delta;
 A(:,i)=mean(seff3(th+d,x,y)-seff3(th-d,x,y),2)/d(i)/2;
end
cov=diag((A\B)/(A')/n);%'

function cov=estcov4_d(th,x,y)
n=length(y);
p=length(th);
f=seff4(th,x,y);
B=f*f'/n;%'  
delta=1e-6;
for i=1:p
 d=zeros(p,1);					
 d(i)=th(i)*delta;
 A(:,i)=mean(seff4(th+d,x,y)-seff4(th-d,x,y),2)/d(i)/2;
end
cov=diag((A\B)/(A')/n);%'

function f=cov_cal(th2,cov2,th0)
lower2=th2-norminv(0.975)*sqrt(diag(cov2));
upper2=th2+norminv(0.975)*sqrt(diag(cov2));
coverage2=zeros(3,1);
for iter=1:1000
    tmp=sign(th0-lower2(:,iter))+sign(upper2(:,iter)-th0);
    for j=1:3
        if(tmp(j)==2)
            coverage2(j)=coverage2(j)+1;
        end
    end
end
f=coverage2;