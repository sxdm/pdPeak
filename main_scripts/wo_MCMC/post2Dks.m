function[ppdf]=post2Dks(muMminusmuF,sigma,lndata)
% calculate 2D posterior pdf of um-uf vs sigma
% assuming flat prior and 1:1 sex ratio.
% data to be a column vector
% muMminusmuF to ba a scaler
% sigma to ba a scaler
k=muMminusmuF;
N=length(lndata);
Sex=dec2bin(0:2^N-1)=='1';
Sex([1,end],:)=[];% The probability of having all M or all F is ignored.
Nm=sum(1-Sex,2);
Nf=sum(Sex,2);
D=lndata'*lndata-((Sex*lndata).^2./Nf+((1-Sex)*lndata).^2./Nm);
Ck=Nm.*Nf.*(k+(Sex*lndata)./Nf-((1-Sex)*lndata)./Nm).^2./N+D;
lpdf=-gammaln(0.5.*N-1.5)+(2-0.5.*N).*log(2)-0.5.*log(pi)-0.5.*log(N)-(N-1).*log(sigma)+...
    log(sum(exp(-0.5.*Ck./(sigma.^2)),1))-log(sum(D.^(1.5-0.5.*N)./sqrt(Nm.*Nf),1));
ppdf=exp(lpdf);
end