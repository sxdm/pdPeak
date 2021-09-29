function[est]=bdi(data)
%Lovejoy et al. 1989, in Hominidae, Jaka Book.
n=length(data);
bino=zeros(n-1,1);
mfs=zeros(n-1,1);
sdata=sort(data,'ascend');
for i=1:(n-1)
    t=sum(log(1:n))-sum(log(1:i))-sum(log(1:(n-i)));
    t=t+log(0.5).*n;
    bino(i)=exp(t);
    mfs(i)=mean(sdata(i+1:end))./mean(sdata(1:i));
end
bino=bino./sum(bino);
est=sum(mfs.*bino);
end