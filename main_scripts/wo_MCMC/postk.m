function[pdf]=postk(muMminusmuF,lndata)
% calculate posterior pdf of um-uf
% assuming flat prior and 1:1 sex ratio.
% data to be a column vector
% muMminusmuF to ba a scaler
k=muMminusmuF;
sz=size(k);
pdf=k(:);
for i=1:prod(sz)
N=length(lndata);
Sex=dec2bin(0:2^N-1)=='1';
Sex([1,end],:)=[];% The probability of having all M or all F is ignored.
Nm=sum(1-Sex,2);
Nf=sum(Sex,2);
D=lndata'*lndata-((Sex*lndata).^2./Nf+((1-Sex)*lndata).^2./Nm);
Ck=Nm.*Nf.*(k(i)+(Sex*lndata)./Nf-((1-Sex)*lndata)./Nm).^2./N+D;
lpdf=-0.5.*(log(pi)+log(N))+gammaln(0.5.*N-1)-gammaln(0.5.*N-1.5)+...
    log(sum(Ck.^(1-0.5.*N),1))-log(sum(D.^(1.5-0.5.*N)./sqrt(Nm.*Nf),1));
pdf(i)=exp(lpdf);
end
pdf=reshape(pdf,sz);
end