function[map]=mMCEM_k(Rk,Rs,data)
mcn=length(Rk);
sRk=sort(Rk,2,'ascend');
uk=sRk(fix(0.95.*mcn));
lk=sRk(fix(0.005.*mcn)+1);
dk=(uk-lk)/50;
edge=lk:dk:uk;
hgm=histcounts(Rk,edge);
[~,idx]=max(hgm);
edge2=0.5.*(edge(1:end-1)+edge(2:end));
map=edge2(idx);
N=length(data);
rec=zeros(1,20);
rec(end)=inf;
rec=[rec(2:end),map];
a=inf;
t=1;
while a>dk/10&&t<100
    Rst=Rs(:,(Rk>=(map-dk))&(Rk<=(map+dk)));
    map=fzero(@(k)Qslope(k,Rst),map);
    rec=[rec(2:end),map];
    a=max(rec)-min(rec);
    t=t+1;
end
function[pdifQ]=Qslope(kt,Rstt)
    m=length(Rstt(1,:));
    Data=repmat(data,1,m);
    Nf=sum(Rstt,1);
    Nm=N-Nf;
    k_hat=pickmean(Data,Rstt==0,1)-pickmean(Data,Rstt==1,1);
    Ddash=N.*sum((data-mean(data)).^2)./(Nf.*Nm)-k_hat.^2;
    pdifQ=sum((kt-k_hat)./((kt-k_hat).^2+Ddash),2);
end
end
function[out]=pickmean(data,tf,dim)
out=sum(data.*tf,dim)./sum(tf,dim);
end