function[est]=mom(data)
%Josephson et al, 1996, Am J Phys Anthropol 100.
u1hat=mean(log(data));
rtu2hat=std(log(data),1);
m4hatz=mean(((log(data)-u1hat)./rtu2hat).^4);
if 1.5-0.5*m4hatz<0
    dhatz=nan;
else
    dhatz=(1.5-0.5.*m4hatz).^0.25;
end
if 1-dhatz^2<0
    shatz=nan;
else
    shatz=sqrt(1-dhatz.^2);
end
dhat=dhatz*rtu2hat;
shat=shatz*rtu2hat;
est=[exp(2.*dhat),sqrt(exp(shat.^2)-1)];
end