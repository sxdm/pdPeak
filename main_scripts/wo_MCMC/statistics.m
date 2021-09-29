clear
close all
load resultMAP.mat

dvk=6;
dvs=6;
ftsz=12;

%
n=1000000;%pseudo-number of MCMCsamples to make color density in the bivariate figure
%comparable with that from MCMC.
mk=mm(data);
ms1=std(data(data>mean(data)))./mean(data(data>mean(data)));
ms2=std(data(data<=mean(data)))./mean(data(data<=mean(data)));
ms=mean([ms1,ms2]);
uk=fzero(@(x1)integral(@(x2)postk(x2,log(data)),0,x1).*2-0.9,mk);
lk=fzero(@(x1)integral(@(x2)postk(x2,log(data)),0,x1).*2-0.1,mk);
us=fzero(@(x1)integral(@(x2)posts(x2,log(data)),0,x1)-0.9,ms);
ls=fzero(@(x1)integral(@(x2)posts(x2,log(data)),0,x1)-0.1,ms);
dk=(uk-lk)/50;
ds=(us-ls)/50;
uk2=lk+(uk-lk)*1.5;
lk2=uk-(uk-lk)*1.5;if lk2<0;lk2=0;end
us2=ls+(us-ls)*1.5;
ls2=us-(us-ls)*1.5;if ls2<0;ls2=0;end
%

%
f1=figure;
plotgramUorL(data,nan(size(data(:,1))),1,f1)
%

%
f31=figure;hold on
area(ku,yku,'Facecolor','b','FaceAlpha',0.2)
plot(ku,yku,'g')
scatter(estmgku,interp1(ku,yku,estmgku),[],'rd','filled')
plot([ku(1),ku(end)],[lvku95,lvku95],'b:')
plot([ku(1),ku(end)],[lvku68,lvku68],'b--')
title('Sexual dimorphism (m/f ratio)')
xlabel('m/f ratio')
ylabel('posterior probability density')
set(gca,'YTickLabel',{});
axis([lk2,uk2,0,1.05.*max(yku)])
digk=floor(log10((exp(uk2)-exp(lk2))/dvk));
tk=ceil(10^(-digk)*exp(lk2))*10^digk;
tkd=round(10^(-digk)*(exp(uk2)-exp(lk2))/dvk)*10^digk;
tk=tk:tkd:exp(uk2);
if digk<0
    ep=['%.',num2str(-digk),'f'];
    tks=compose(ep,tk);
else
    tks=compose('%g',tk./10.^digk);
    if digk>0
    text(uk2,0,['x10^',num2str(digk)],'FontSize',ftsz)
    end
end
set(gca,'Xtick',log(tk),'XGrid','on','XTickLabel',tks)
hold off
%

%
f4=figure;hold on
bxul=0.15;
set(gca,'FontSize',ftsz)
area(s,ys,'Facecolor','b','FaceAlpha',0.2)
plot(s,ys,'g')
scatter(estmgs,interp1(s,ys,estmgs),[],'rd','filled')
plot([s(1),s(end)],[lvs95,lvs95],'b:')
plot([s(1),s(end)],[lvs68,lvs68],'b--')
title('Within-sex CV')
xlabel('within-sex CV [%]')
ylabel('posterior probability density')
set(gca,'YTickLabel',{});
axis([ls2,us2,0,1.05.*max(ys)])
digs=floor(log10((sqrt(exp(us2^2)-1)-sqrt(exp(ls2^2)-1))/dvs));
ts=ceil(10^(-digs)*sqrt(exp(ls2^2)-1))*10^digs;
tsd=round(10^(-digs)*(sqrt(exp(us2^2)-1)-sqrt(exp(ls2^2)-1))/dvs)*10^digs;
ts=ts:tsd:sqrt(exp(us2^2)-1);
if digs+2<0
    if digs+2<-1
    t2digs=floor(log10(ts(2)));
    ep=['%.',num2str(t2digs-digs),'f'];
    tss=compose(ep,ts./10.^t2digs);
    if t2digs+2<0
    text(us2,0,['x10^-^',num2str(-(t2digs+2))],'FontSize',ftsz)
    else
    text(us2,0,['x10^',num2str(t2digs+2)],'FontSize',ftsz)
    end
    else
    ep=['%.',num2str(-(digs+2)),'f'];
    tss=compose(ep,ts.*100);
    end
else
    tss=compose('%g',ts./10.^digs);
    if digs+2>0
    text(us2,0,['x10^',num2str(digs+2)],'FontSize',ftsz)
    end
end
set(gca,'Xtick',sqrt(log(ts.^2+1)),'XGrid','on','XTickLabel',tss)
hold off
%

%
f51=figure;
frmw=ls2:ds:us2+ds;
frmr=lk2:dk:uk2+dk;
ss1=0.5.*(frmw(1:end-1)+frmw(2:end));
kk1=0.5.*(frmr(1:end-1)+frmr(2:end));
[X2,Y2]=meshgrid(ss1,kk1);
frame=X2;
frame(:,:)=0;
for ii=1:length(X2(1,:))
    for jj=1:length(Y2(:,1))
        frame(jj,ii)=post2Dks(Y2(jj,ii),X2(jj,ii),log(data)).*2.*dk.*ds.*n;
    end
end
set(gca,'FontSize',ftsz)
imagesc(frame,'XData',frmw,'YData',frmr)
col=[(1:-1/63:0)',(1:-1/63:0)',ones(64,1)];
colormap(col)
caxis([0,800])
set(gca,'YDir','normal')
if digs+2<-1
    if t2digs+2<0
    text(frmw(end),frmr(1),['x10^-^',num2str(-(t2digs+2))],'FontSize',ftsz)
    else
    text(frmw(end),frmr(1),['x10^',num2str(t2digs+2)],'FontSize',ftsz)
    end
elseif digs+2>0
    text(frmw(end),frmr(1),['x10^',num2str(digs+2)],'FontSize',ftsz)
end
set(gca,'Xtick',sqrt(log(ts.^2+1)),'XTickLabel',tss)
set(gca,'Ytick',log(tk),'YTickLabel',tk)
axis square
hold on
set(gca,'FontSize',ftsz)
[X1,Y1]=meshgrid(frmw,frmr);
CV2d=sqrt(2.*(1+exp(Y1).^2).*(sqrt(exp(X1.^2)-1).^2+1)./(1+exp(Y1)).^2-1).*100;
cinc=sqrt(2.*(1+exp(frmr(1)).^2).*(ts.^2+1)./(1+exp(frmr(1))).^2-1).*100;
cinc=[cinc,cinc(end)+mean(diff(cinc)):mean(diff(cinc)):max(CV2d(:))];
if length(cinc)>30
    cinc=min(CV2d(:)):(max(CV2d(:))-min(CV2d(:)))/5:max(CV2d(:));
end
[Ct2,h2]=contour(X1,Y1,CV2d,cinc,'k:');
set(h2,'ShowText','on')
clabel(Ct2,h2,'FontSize',0.8.*ftsz)
title({'Bivariate dinsity plot','(dotted contour = total cv level)'})
xlabel('within-sex cv [%]')
ylabel('m/f ratio')
axis([frmw(1),frmw(end),frmr(1),frmr(end)])
hold off
axis square

%
ci68cv=cell(length(l68s),1);
for i=1:length(l68s)
    ci68cv{i}=[num2str(sqrt(exp(l68s(i)^2)-1),'%10.3f'),'-',num2str(sqrt(exp(u68s(i)^2)-1),'%10.3f'),' '];
end
t=ci68cv{1};
for i=2:length(l68s)
    t=[t,' ',ci68cv{i}];
end
ci68cv=t;

ci95cv=cell(length(l95s),1);
for i=1:length(l95s)
    ci95cv{i}=[num2str(sqrt(exp(l95s(i)^2)-1),'%10.3f'),'-',num2str(sqrt(exp(u95s(i)^2)-1),'%10.3f'),' '];
end
t=ci95cv{1};
for i=2:length(l95s)
    t=[t,' ',ci95cv{i}];
end
ci95cv=t;

ci68sdmu=cell(length(l68ku),1);
for i=1:length(l68ku)
    ci68sdmu{i}=[num2str(exp(l68ku(i)),'%10.3f'),'-',num2str(exp(u68ku(i)),'%10.3f'),' '];
end
t=ci68sdmu{1};
for i=2:length(l68ku)
    t=[t,' ',ci68sdmu{i}];
end
ci68sdmu=t;

ci95sdmu=cell(length(l95ku),1);
for i=1:length(l95ku)
    ci95sdmu{i}=[num2str(exp(l95ku(i)),'%10.3f'),'-',num2str(exp(u95ku(i)),'%10.3f'),' '];
end
t=ci95sdmu{1};
for i=2:length(l95ku)
    t=[t,' ',ci95sdmu{i}];
end
ci95sdmu=t;

subjects={' ';...
    'm/f ratio (by pdPeak)   = ';...
    '  68% credible interval = ';...
    '  95% credible interval = ';...
    'w-sx CV   (by pdPeak)   = ';...
    '  68% credible interval = ';...
    '  95% credible interval = ';...
    'by   mean method        = ';...
    'by   BDI  method        = ';...
    'by   MOM  method        = ';...
    'by   CVM(Plavcan1994)   = ';...
    ' ';...
    'Sample statistics';...
    '                  N = ';...
    '               mean = ';...
    '                 sd = ';...
    '             totCV  = ';...
    '             totCV* = ';...
    'skewness (log-base) = ';...
    '       p-val (skew) = ';...
    ' '};
values={' ';...
    num2str(exp(estmgku),'%10.3f');...
    ci68sdmu;...
    ci95sdmu;...
    num2str(sqrt(exp(estmgs.^2)-1),'%10.3f');...
    ci68cv;...
    ci95cv;...
    num2str(mm(data(~isnan(data))),'%10.3f');...
    num2str(bdi(data(~isnan(data))),'%10.3f');...
    num2str(mom(data(~isnan(data))),'%10.3f');...
    num2str(lcv(data(~isnan(data))),'%10.3f');...
    '';...
    '';...
    num2str(length(data));...
    num2str(mean(data));...
    num2str(std(data,0));...
    num2str(std(data,0)./mean(data));...
    [num2str(std(data,0)./mean(data).*(1+0.25./length(data))),' *Sokal&Braumann(1980) correction'];...
    num2str(skw1(log(data)));...
    [num2str(DAgostino(log(data),'both')),' (two-tailed)'];...
    ' '};

for i=1:length(subjects)
    disp([subjects{i},values{i}])
end