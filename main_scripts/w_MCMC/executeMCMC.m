clear
close all

load tmp.mat
data=data(data>0);
datal=log(data);

ini=10000;%length of burn-in period
t=1010000;%chain length
imuU=mean(datal(~isnan(datal(:,1)),1));
init_mu=[mean(datal(datal(:,1)>=imuU,1));mean(datal(datal(:,1)<imuU,1))];
init_s=repmat(mean([std(datal(datal(:,1)>=imuU,1));std(datal(datal(:,1)<imuU,1))]),2,1);
init_logicals=datal(:,1)<imuU;

[Qm,Qs,S]=sxdmMCMC5(t,init_mu,init_s,init_logicals,datal);

Rm=Qm(:,ini+1:end);
RS=S(:,ini+1:end);
Rs=Qs(:,ini+1:end);

clearvars -except Rm RS Rs data ini t

save('resultMCMC.mat')