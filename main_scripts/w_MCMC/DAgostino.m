function[pval,varargout]=DAgostino(data,oneorboth)
%D'Agostino test: D'Agostino (1970) Biometrika 57(3):679-681.
data=data(:);
n=length(data);
if n<8
    pval=nan;
    varargout{1}=nan;
else
xbar=mean(data);
g1=(sum((data-xbar).^3)/n)/sqrt(sum((data-xbar).^2)/n)^3;
Y=g1*sqrt((n+1)*(n+3)/6/(n-2));
b2=3*(n^2+27*n-70)*(n+1)*(n+3)/(n-2)/(n+5)/(n+7)/(n+9);
W=sqrt(-1+sqrt(2*(b2-1)));
s=1/sqrt(log(W));
a=sqrt(2/(W^2-1));
z=s*log(Y/a+sqrt((Y/a)^2+1));
if z>0
    pval=1-0.5.*(1+erf(z/sqrt(2)));
elseif z<=0
    pval=0.5.*(1+erf(z/sqrt(2)));
end
varargout{1}=z;
end
if strcmp(oneorboth,'both')
    pval=2.*pval;
end
if isnan(pval)
    pval='insufficient N';
end
end