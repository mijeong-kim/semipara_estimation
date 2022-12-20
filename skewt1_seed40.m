function skewt1_seed40() 

path='/usr/local/bin/R';
libpath='/Library/Frameworks/R.framework/Resources/library';
Rinit('sn',path,libpath)

n=300; simus=1000; beta0=[5, 1, 1.8];
xi=0; omega=3; alpha=3; nu=10; p=0.7;
bv=sqrt(nu)*gamma(0.5*(nu-1))/sqrt(pi)/gamma(0.5*nu);
d=alpha/sqrt(1+alpha^2);
%mu1_1=xi+omega*bv*d;
%Gamma(a,b)
a=2.5;b=3;[a,b]
mu2_1=omega^2*(nu/(nu-2)-(bv*d)^2);
mu2_2=a*b^2;
mu2=p*mu2_1+(1-p)*mu2_2;

theta0=[beta0, mu2];

Rrun('load("skt1_seed40.RData")')
[xy,b_olse,result1,coverage1]=Rpull('xy','b_olse','result','coverage');


th4=zeros(4,simus); cov4=zeros(4,simus); 
th4_eflag=zeros(1,simus);  
samplevar=zeros(2,simus);

%options2=optimset('MaxFunEvals',1e4,'display','off');
options1=optimset('algorithm','levenberg-marquardt','MaxFunEvals',1e6,'MaxIter',1e6,'display','off');
options3=optimset('algorithm','levenberg-marquardt','MaxFunEvals',1e4,'MaxIter',1e6,'display','off');
%options4=optimset('algorithm','trust-region-reflective','MaxFunEvals',1e4,'MaxIter',1e6,'display','off');
%options3=optimset('MaxFunEvals',1e4,'MaxIter',1e4,'display','off');
%options6 = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1500)


fails4=0;  

for iter=1:simus
 	iter
    eflag=1; 
     x=xy(:,1:2,iter)';
     y=xy(:,3,iter)';
    %STARTING Semiparametric estimator
    if (eflag>0)
     [th4(:,iter),f,eflag]=fsolve(@seffscaled4,b_olse(iter,:)',options3,x,y);
     cov4(:,iter)=estcov4_d(th4(:,iter),x,y);
     e4=y-m(x,th4(:,iter));
     th4_eflag(iter)=eflag;
     samplevar(1,iter)=var(e4);
     if (eflag<=0)
       eflag=0;
       fails4=fails4+1;
     end
   end
end

fails4

'seed40'

coverage1 
result1

coverage4=cov_cal(th4,cov4,theta0,simus)/simus*100;
coverage4'
result2=[median(th4,2) std(th4')'   median(sqrt(cov4),2)]



save skewt1_seed40.mat result1 coverage1  result2 th4  cov4  samplevar fails4 coverage4 th4_eflag 
Rclear;



function f=e2sum(th,x,y)
f=sum((y-m(x,th)).^2);


function f=m(x,beta)
f=beta(1)+beta(2)*x(1,:)+beta(3)*x(2,:);


function f=dm(x,beta)
f=[ones(1,length(x));x];


function f=seff4(th,x,y)
e=y-m(x,th); 
n=size(e,2);
%h=0.9*min(std(e),iqr(e)/1.34)*n^(-1/5);
h=1.06*std(e)*n^(-1/5);
mu2=th(4);
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




function f=seffscaled4(th,x,y)
f=mean(seff4(th,x,y),2);

function cov=estcov4(th,x,y)
n=length(y);
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


function f=cov_cal(th2,cov2,th0,simus)
lower2=th2-norminv(0.975)*sqrt(diag(cov2));
upper2=th2+norminv(0.975)*sqrt(diag(cov2));
p=size(th2,1);
coverage2=zeros(p,1);
for iter=1:simus
    tmp=sign(th0'-lower2(:,iter))+sign(upper2(:,iter)-th0');
    for j=1:p
        if(tmp(j)==2)
            coverage2(j)=coverage2(j)+1;
        end
    end
end
f=coverage2;
