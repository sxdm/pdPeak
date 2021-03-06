function[map]=mMCEM_s(RS,Rs,data)
mcn=length(RS);
sRS=sort(RS,2,'ascend');
us=sRS(fix(0.95.*mcn));
ls=sRS(fix(0.005.*mcn)+1);
ds=(us-ls)/50;
edge=ls:ds:us;
hgm=histcounts(RS,edge);
[~,idx]=max(hgm);
edge2=0.5.*(edge(1:end-1)+edge(2:end));
map=edge2(idx);
N=length(data);
rec=zeros(1,20);
rec(end)=inf;
rec=[rec(2:end),map];
a=inf;
t=1;
while a>ds/10&&t<100
    Rst=Rs(:,(RS>=(map-ds))&(RS<=(map+ds)));
    map=mxQs(Rst);
    rec=[rec(2:end),map];
    a=max(rec)-min(rec);
    t=t+1;
end
function[st]=mxQs(Rstt)
    m=length(Rstt(1,:));
    Data=repmat(data,1,m);
    Nf=sum(Rstt,1);
    Nm=N-Nf;
    D=sum(data.^2)-(((sum(Data.*(Rstt==0),1)).^2)./Nm+((sum(Data.*(Rstt==1),1)).^2)./Nf);
    st=sqrt(sum(D,2)/(m*(N-2)));
end
end