function[intvls,varargout]=findxxintvlk_MC(Rk,percent,smoothitr)
%
frck=30;
%
Rk=abs(Rk);
n=length(Rk);
MCsample=sort(Rk,2,'ascend');
u5=MCsample(fix(n.*0.95));
l=MCsample(fix(n.*0.005)+1);
dk=(u5-l)/frck;
mx=max(Rk);
k=0:dk:(mx+dk);
if k(end-1)~=mx
    k(end)=[];
end
fq=histcounts(Rk,k);
pdf=fq./(n.*diff(k));
pdf=[pdf(end:-1:1),pdf];
k2=0.5.*(k(1:end-1)+k(2:end));
k2=[-k2(end:-1:1),k2];

if smoothitr~=0
    k3=[k2(1),0.5.*(k2(1:end-1)+k2(2:end)),k2(end)];
    for i=1:smoothitr
        pdf3=interp1(k2,pdf,k3);
        pdf=interp1(k3,pdf3,k2);
    end
end

pdf0=interp1(k2,pdf,0);
pdf=[pdf0,pdf(k2>0)];
k2=[0,k2(k2>0)];
[intvls,lv]=findxxintvl(k2,pdf,percent,dk/1000);
varargout{1}=[k2;pdf];
varargout{2}=lv;
end
function[xxintvls,y]=findxxintvl(x,pdf,percent,dx)
pdf=pdf./calarea([x(1),x(end)]);
[mapd,~]=max(pdf);
pfun=@(lv)calarea(findintvl(lv));
y=fzero(@(xxx)pfun(xxx)-percent/100,[mapd,0]);
xxintvls=findintvl(y);
    function[area]=calarea(intvls)
        if isempty(intvls)
            area=0;
        else
            ni=length(intvls(:,1));
            areas=zeros(ni,1);
            for j=1:ni
                areas(j)=Integ(@(xx)interp1(x,pdf,xx),intvls(j,1),intvls(j,2),dx);
            end
            area=sum(areas);
        end
    end
    function[intvls]=findintvl(level)
        t=pdf>=level;
        d=diff(t);
        if t(1)&&t(end)
            n=sum(d==1)+1;
            intvls=nan(n,2);
            idx1=find(d==1);
            idx_1=find(d==-1);
            if n==1
                intvls=[x(1),x(end)];
            else
                intvls(1,1)=[x(1),interp1([pdf(idx_1(1)),pdf(idx_1(1)+1)],[x(idx_1(1)),x(idx_1(1)+1)],level)];
                intvls(end,2)=[interp1([pdf(idx1(end)),pdf(idx1(end)+1)],[x(idx1(end)),x(idx1(end)+1)],level),x(end)];
                if n>2
                    for i=2:(n-1)
                        intvls(i,:)=[interp1([pdf(idx1(i-1)),pdf(idx1(i-1)+1)],[x(idx1(i-1)),x(idx1(i-1)+1)],level),...
                            interp1([pdf(idx_1(i)),pdf(idx_1(i)+1)],[x(idx_1(i)),x(idx_1(i)+1)],level)];
                    end
                end
            end
        elseif t(1)
            n=sum(d==1)+1;
            intvls=nan(n,2);
            idx1=find(d==1);
            idx_1=find(d==-1);
            intvls(1,:)=[x(1),interp1([pdf(idx_1(1)),pdf(idx_1(1)+1)],[x(idx_1(1)),x(idx_1(1)+1)],level)];
            if n>1
                for i=2:n
                    intvls(i,:)=[interp1([pdf(idx1(i-1)),pdf(idx1(i-1)+1)],[x(idx1(i-1)),x(idx1(i-1)+1)],level),...
                        interp1([pdf(idx_1(i)),pdf(idx_1(i)+1)],[x(idx_1(i)),x(idx_1(i)+1)],level)];
                end
            end
        elseif t(end)
            n=sum(d==1);
            intvls=nan(n,2);
            idx1=find(d==1);
            idx_1=find(d==-1);
            if n>1
                for i=1:(n-1)
                    intvls(i,:)=[interp1([pdf(idx1(i)),pdf(idx1(i)+1)],[x(idx1(i)),x(idx1(i)+1)],level),...
                        interp1([pdf(idx_1(i)),pdf(idx_1(i)+1)],[x(idx_1(i)),x(idx_1(i)+1)],level)];
                end
            end
            intvls(end,:)=[interp1([pdf(idx1(end)),pdf(idx1(end)+1)],[x(idx1(end)),x(idx1(end)+1)],level),x(end)];
        else
            n=sum(d==1);
            intvls=nan(n,2);
            idx1=find(d==1);
            idx_1=find(d==-1);
            if n>0
                for i=1:n
                    intvls(i,:)=[interp1([pdf(idx1(i)),pdf(idx1(i)+1)],[x(idx1(i)),x(idx1(i)+1)],level),...
                        interp1([pdf(idx_1(i)),pdf(idx_1(i)+1)],[x(idx_1(i)),x(idx_1(i)+1)],level)];
                end
            end
        end
    end
end
function [ S ] = Integ( func, a, b, dx )
f=fcnchk(func);
sep=[a:dx:b,b];
dxs=diff(sep);
input=0.5.*(sep(1:end-1)+sep(2:end));
output=f(input);
S=sum(output.*dxs);
end