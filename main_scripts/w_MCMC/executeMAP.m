clear
close all

load resultMCMC.mat

%
n=length(Rm(1,:));
%

%
Rku=Rm(1,:)-Rm(2,:);
Rst=Rs;
Rst(:,Rku<0)=(1-Rst(:,Rku<0));
Rkut=abs(Rku);
estmgku=mMCEM_k(Rkut,Rst,log(data));
estmgs=mMCEM_s(RS(1,:),Rs,log(data));
%

%
[ci95,ys,lvs95]=findxxintvls_MC(RS(1,:),95,3);
[ci68,~,lvs68]=findxxintvls_MC(RS(1,:),68,3);
u95s=ci95(:,2);
l95s=ci95(:,1);
u68s=ci68(:,2);
l68s=ci68(:,1);
s=ys(1,:);
ys=ys(2,:);

[ci95,yku,lvku95]=findxxintvlk_MC(Rku,95,3);
[ci68,~,lvku68]=findxxintvlk_MC(Rku,68,3);
u95ku=ci95(:,2);
l95ku=ci95(:,1);
u68ku=ci68(:,2);
l68ku=ci68(:,1);
ku=yku(1,:);
yku=yku(2,:);

%

clearvars -except estmgku estmgs s ku yku ys...
    l68s u68s l95s u95s l68ku u68ku l95ku u95ku...
    lvku95 lvku68 lvs95 lvs68

save('resultMAP.mat')