clear
close all

load tmp.mat
data=data(data>0);

%
estmg=abs(map_marg(log(data)));
estmgku=estmg(1);
estmgs=estmg(2);
%

%
[ci95,ys,lvs95]=findxxintvl2s(95,log(data));
[ci68,~,lvs68]=findxxintvl2s(68,log(data));
u95s=ci95(:,2);
l95s=ci95(:,1);
u68s=ci68(:,2);
l68s=ci68(:,1);
s=ys(1,:);
ys=ys(2,:);

[ci95,yku,lvku95]=findxxintvl2k(95,log(data));
[ci68,~,lvku68]=findxxintvl2k(68,log(data));
u95ku=ci95(:,2);
l95ku=ci95(:,1);
u68ku=ci68(:,2);
l68ku=ci68(:,1);
ku=yku(1,:);
yku=yku(2,:);
%

clearvars -except estmgku estmgs s ku yku ys...
    l68s u68s l95s u95s l68ku u68ku l95ku u95ku...
    lvku95 lvku68 lvs95 lvs68...
    data

save('resultMAP.mat')