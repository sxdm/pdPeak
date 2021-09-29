function[out]=skw1(data)
data(isnan(data))=[];
mu=mean(data);
sd=std(data,1);
out=mean((data-mu).^3)./sd.^3;
end