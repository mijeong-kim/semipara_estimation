function table3_n300() %Gaussian mixture
for i=12
randn('state',i); rand('state',i);
simus=1000;n=300; 
mu=[-2,3]; sigma=[0.6, 0.7]; 
p1=0.6; p2=1-p1; 
m1=p1*mu(1)+p2*mu(2);
mu2=p1*(mu(1)^2+sigma(1)^2)+p2*(mu(2)^2+sigma(2)^2)-m1^2;
beta0=[12 -0.5]';
theta0=[beta0; mu2];
th1=zeros(2,simus); 
th2=zeros(3,simus); th3=zeros(3,simus); th4=zeros(3,simus);
cov2=zeros(3,simus); cov3=zeros(3,simus); cov4=zeros(3,simus);
samplevar=zeros(4,simus);
e1_set=zeros(n,simus);  e2_set=zeros(n,simus);  e3_set=zeros(n,simus);  e4_set=zeros(n,simus);
H=zeros(1,simus); pValue=zeros(1,simus); W=zeros(1,simus);

%options2=optimset('MaxFunEvals',1e4,'display','off');
options1=optimset('algorithm','levenberg-marquardt','MaxFunEvals',1e6,'MaxIter',1e6,'display','off');
options3=optimset('algorithm','levenberg-marquardt','MaxFunEvals',1e4,'MaxIter',1e6,'display','off');
%options4=optimset('algorithm','trust-region-reflective','MaxFunEvals',1e4,'MaxIter',1e6,'display','off');
%options3=optimset('MaxFunEvals',1e4,'MaxIter',1e4,'display','off');
options5 = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1e6);
%options6 = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1500)
i
iter=1;
fails1=0;fails2=0; fails3=0; fails4=0;
while (iter<=simus)
  [x,y]=gendata(mu,sigma,p1,beta0,n);
 eflag=1;
 %OLS
 if (eflag>0)
   [th1(:,iter),f,eflag]=fminsearch(@e2sum,beta0,[],x,y);
    e1=y-m(x,th1(:,iter));
    if (eflag<=0) 
     [iter 1 eflag th1(:,iter)']%'
     fails1=fails1+1;
   end
 end
%STARTING Semiparametric estimator
 if (eflag>0)
     [th2(:,iter),f,eflag]=fsolve(@seffscaled,theta0,options1,x,y);
     cov2(:,iter)=estcov_d(th2(:,iter),x,y);
     e2=y-m(x,th2(:,iter));
     if (eflag<=0)
       eflag=0;
       fails2=fails2+1;
     end
 end
  if (eflag>0)
     [th3(:,iter),f,eflag]=fsolve(@seffscaled3,theta0,options3,x,y);
     cov3(:,iter)=estcov3_d(th3(:,iter),x,y);
     e3=y-m(x,th3(:,iter));
     if (eflag<=0)'
       eflag=0;
       fails3=fails3+1;
     end
  end
    if (eflag>0)
     [th4(:,iter),f,eflag]=fsolve(@seffscaled4,theta0,options3,x,y);
     cov4(:,iter)=estcov4_d(th4(:,iter),x,y);
     e4=y-m(x,th4(:,iter));
     if (eflag<=0)
       eflag=0;
       fails4=fails4+1;
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
[fails1 fails2 fails3 fails4]


coverage2=cov_cal(th2,cov2,theta0)/simus*100
coverage3=cov_cal(th3,cov3,theta0)/simus*100
coverage4=cov_cal(th4,cov4,theta0)/simus*100

[median(th2,2) std(th2')'   median(sqrt(cov2),2)]
[median(th3,2) std(th3')'   median(sqrt(cov3),2)]
[median(th4,2) std(th4')'   median(sqrt(cov4),2)]
save table3_n300.mat th1 th2 th3 th4 cov2 cov3 cov4 samplevar i fails1 fails2 fails3 fails4 e1_set e2_set e3_set e4_set H pValue W n coverage2 coverage3 coverage4 
end


function [x,y]=gendata(mu,sigma,p,beta,n)
%x=unifrnd(0,7,1,n); 
%x=chi2rnd(4.5,1,n);
pd1= makedist('Gamma','a',2.5,'b',1.5);
x=random(pd1,1,n);
e=zeros(1,n);
m1=mu(1); m2=mu(2); sigma1=sigma(1); sigma2=sigma(2);
for i=1:n
 u=unifrnd(0,1,1);
 if (u<p)
     e(i)=normrnd(m1,sigma1,1);
 else
     e(i)=normrnd(m2,sigma2,1);
 end
end   
m0=p*mu(1)+(1-p)*mu(2);
e1=e-m0;
y=m(x,beta)+e1;

function f=e2sum(th,x,y)
f=sum((y-m(x,th)).^2);

function f=f1f(e,mu,sigma,p)
m1=mu(1); m2=mu(2); 
sigma1=sigma(1); sigma2=sigma(2);
p1=p(1); p2=p(2);
t1=(e-m1)/sigma1; t2=(e-m2)/sigma2;   
f1=1/sigma1*exp(-t1.^2/2);
f2=1/sigma2*exp(-t2.^2/2);
B=-p1*t1.*f1/sigma1-p2*t2.*f2/sigma2;
A=p1*f1+p2*f2;
f=B./A;

function y=m(x,beta)
y=beta(1)*exp(beta(2)*x);

function y=dm(x,theta)
y=[exp(theta(2)*x) 
    m(x,theta).*x];


function f=seff(th,x,y)
options = statset('MaxIter',1e6);
%options3=optimset('algorithm','levenberg-marquardt','MaxFunEvals',1e4,'MaxIter',1e6,'display','off');
e=y-m(x,th);
mu2=th(3);
result=fitgmdist(e',2,'CovarianceType','diagonal','Options',options,'RegularizationValue',0.1);
p=result.PComponents;
mu=result.mu;
sigma=sqrt(result.Sigma);
mu3=mu_3(mu,sigma,p);
mu4=mu_4(mu,sigma,p);
t=e.^2-mu2-mu3*e./mu2;
Et2=mu4-mu2^2-mu3^2/mu2;
f=[-f1f(e,mu,sigma,p).*(dm(x,th)-mean(dm(x,th),2))+(e./mu2-mu3*t./(mu2*Et2)).*mean(dm(x,th),2);
    t./Et2];

function f=mu_3(mu,sigma,p)
m1=mu(1); m2=mu(2); 
sigma1=sigma(1); sigma2=sigma(2);
mu3_1=3*m1*sigma1^2+m1^3;
mu3_2=3*m2*sigma2^2+m2^3;
f=p(1)*(mu3_1)+p(2)*mu3_2;

function f=mu_4(mu,sigma,p)
m1=mu(1); m2=mu(2); 
sigma1=sigma(1); sigma2=sigma(2);
mu4_1=3*sigma1^4+6*m1^2*sigma1^2+m1^4;
mu4_2=3*sigma2^4+6*m2^2*sigma2^2+m2^4;
f=p(1)*(mu4_1)+p(2)*mu4_2;


function f=seff3(th,x,y)
e=y-m(x,th);  
mu2=th(3);
Et2=2*mu2^2;
f=[e./mu2.*dm(x,th);
    (e.^2-mu2)./Et2];

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
%mu2=sum(e.^2+h^2);
%mu3=sum(e.^3+3*e*h^2)/n;
mu3=sum(e.^3)/n;
%mu4=sum(e.^4+6*e.^2*h^2+3*h^4)/n;
mu4=sum(e.^4)/n;
t=e.^2-mu2-mu3.*e./mu2;
Et2=mu4-mu2.^2-mu3.^2./mu2;
%t2=mean((e.^2-mu2).^2);
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

function cov=estcov(th,x,y)
n=length(y);
p=length(th);
f=seff(th,x,y);
B=f*f'/n;%' 
cov=diag(inv(B)/n);

function cov=estcov3_d(th,x,y)
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

function cov=estcov4(th,x,y)
n=length(y);
p=length(th);
f=seff4(th,x,y);
B=f*f'/n;%' 
cov=diag(inv(B)/n);

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
