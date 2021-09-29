function[]=MCcheck(Qm,S)
t=length(Qm(1,:));
set(0,'Units','pixels') 
scnsize=get(0,'ScreenSize');
pos1=[scnsize(3)/32,scnsize(4)*2/64,scnsize(3)*14/32,scnsize(4)*60/64];
pos2=[scnsize(3)*16/32,scnsize(4)*2/64,scnsize(3)*14/32,scnsize(4)*60/64];

f1=figure;
set(f1,'OuterPosition',pos1)
subplot(3,1,1)
plot(1:t,Qm(1,:))
title({'chain profiles',' '})
xlabel('sample# in chain (burn-in period excluded)')
ylabel('mu for male (muM)')
subplot(3,1,2)
plot(1:t,Qm(2,:))
xlabel('sample# in chain (burn-in period excluded)')
ylabel('mu for female (muF)')
subplot(3,1,3)
plot(1:t,S(1,:))
xlabel('sample# in chain (burn-in period excluded)')
ylabel('sigma')
autoCn=100;
crfQm=zeros(2,autoCn);
crfS=zeros(1,autoCn);
for i=1:autoCn
    tt=corrcoef(Qm(1,1:end-i),Qm(1,i+1:end));
    crfQm(1,i)=tt(2,1);
    tt=corrcoef(Qm(2,1:end-i),Qm(2,i+1:end));
    crfQm(2,i)=tt(2,1);
    tt=corrcoef(S(1,1:end-i),S(1,i+1:end));
    crfS(1,i)=tt(2,1);
end
f2=figure;
set(f2,'OuterPosition',pos2)
subplot(3,1,1)
plot(1:autoCn,crfQm(1,:))
title({'auto-correlation','=correlation coefficient btwn pairs of samples with a distance within the chain'})
xlabel('distance')
ylabel('auto-correlation (muM)')
axis([1,autoCn,-0.3,0.3])
subplot(3,1,2)
plot(1:autoCn,crfQm(2,:))
xlabel('distance')
ylabel('auto-correlation (muF)')
axis([1,autoCn,-0.3,0.3])
subplot(3,1,3)
plot(1:autoCn,crfS(1,:))
xlabel('distance')
ylabel('auto-correlation (sigma)')
axis([1,autoCn,-0.3,0.3])
end