function [Qm,Qs,S] = sxdmMCMC5(t,init_mu,init_s,init_logicals,data_ok)
Qm=zeros(length(init_mu),t);
Qs=zeros(length(init_logicals),t);
S=zeros(2,t);
Qm(:,1)=init_mu;
Qs(:,1)=init_logicals;
S([1,2],1)=init_s;
n=length(init_logicals);
alpha=0.5.*(n-1);

for i=1:t-1
dataad=data_ok;


%Gibbs sampling for the means
dataadM=dataad(Qs(:,i)==0,:);
dataadF=dataad(Qs(:,i)==1,:);
Mn=sum(Qs(:,i)==0);
Fn=sum(Qs(:,i)==1);
mm=mean(dataadM,1);
fm=mean(dataadF,1);
Qm(1,i+1)=randn(1).*S(1,i)./sqrt(Mn)+mm;
Qm(2,i+1)=randn(1).*S(2,i)./sqrt(Fn)+fm;


%Gibbs sampling for the sigma
mvec=dataadM-Qm(1,i+1);
fvec=dataadF-Qm(2,i+1);
ms=mvec'*mvec;
fs=fvec'*fvec;
beta=0.5.*(ms+fs);
S([1,2],i+1)=sqrt(rndinvg(1,alpha,beta));


%Gibbs sampling from multinomial distribution for sex assignments
cpostpM=exp(-0.5.*((dataad-Qm(1,i+1))./S(1,i+1)).^2)./S(1,i+1);
cpostpF=exp(-0.5.*((dataad-Qm(2,i+1))./S(2,i+1)).^2)./S(2,i+1);
cpostpM=cpostpM./(cpostpM+cpostpF);
k=0;
while k==0
    Qstemp=rand(n,1)>cpostpM;
    if ~(all(Qstemp)||all(Qstemp==0))
        k=1;
    end
end
Qs(:,i+1)=Qstemp;


%Random permutation
pyn=rand>0.5;
if pyn
    Qm(:,i+1)=Qm([2,1],i+1);
    S(:,i+1)=S([2,1],i+1);
    Qs(:,i+1)=1-Qs(:,i+1);
end


%figuring
disp([num2str(i),'/',num2str(t)])
end

end
function[output]=rndinvg(n,alpha,beta)
uni=rand(n);
output=1./gammaincinv(uni,alpha);
tmp=1./(gammaincinv(1-uni,alpha,'upper'));
output(uni<0.1)=tmp(uni<0.1);
output=beta.*output;
end