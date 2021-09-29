function[pdf]=posts(sigma,lndata)
% calculate posterior pdf of um-uf
% assuming flat prior and 1:1 sex ratio.
% data to be a column vector
% sigma to ba a scaler
sz=size(sigma);
pdf=sigma(:);
for i=1:prod(sz)
if sigma(i)<=0
    pdf(i)=0;
else
N=length(lndata);
Sex=dec2bin(0:2^N-1)=='1';
Sex([1,end],:)=[];% The probability of having all M or all F is ignored.
Nm=sum(1-Sex,2);
Nf=sum(Sex,2);
D=lndata'*lndata-((Sex*lndata).^2./Nf+((1-Sex)*lndata).^2./Nm);
lpdf=(2.5-0.5.*N).*log(2)-gammaln(0.5.*N-1.5)-(N-2).*log(sigma(i))+...
    log(sum(exp(-0.5.*D./sigma(i).^2)./sqrt(Nm.*Nf),1))-log(sum(D.^(1.5-0.5.*N)./sqrt(Nm.*Nf),1));
pdf(i)=exp(lpdf);
end
end
pdf=reshape(pdf,sz);
end