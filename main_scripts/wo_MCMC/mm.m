function[est]=mm(data)
pvt=mean(data);
hm=data(data>pvt);
hf=data(data<=pvt);
est=mean(hm)./mean(hf);
end