clear
close all
load resultMCMC.mat
load resultMAP.mat

dvk=6;
dvs=6;
ftsz=12;

%
n=length(Rm(1,:));
Rku=Rm(1,:)-Rm(2,:);
Rku1=sort(abs(Rku),'ascend');
uk=Rku1(ceil(0.9.*n));
lk=Rku1(ceil(0.1.*n));
RS1=sort(RS(1,:),'ascend');
us=RS1(ceil(0.9.*n));
ls=RS1(ceil(0.1.*n));
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
MCcheck(Rm,RS)
%

%
f31=figure;hold on
dbox=dk;
box=lk2-dk:dk:uk2+dk;
Rkut=abs(Rku);
B=histcounts(Rkut,box);
set(gca,'FontSize',ftsz)
bar(0.5.*(box(2:end)+box(1:(end-1))),B)
plot(ku,yku.*n.*dbox,'g')
scatter(estmgku,interp1(ku,yku.*n.*dbox,estmgku),[],'rd','filled')
plot([ku(1),ku(end)],[lvku95.*n.*dbox,lvku95.*n.*dbox],'b:')
plot([ku(1),ku(end)],[lvku68.*n.*dbox,lvku68.*n.*dbox],'b--')
title('Sexual dimorphism (m/f ratio)')
xlabel('m/f ratio')
axis([lk2,uk2,0,1.05.*max(B)])
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
dbox=ds;
box=ls2-ds:ds:us2+ds;
B=histcounts(RS(1,:),box);
set(gca,'FontSize',ftsz)
bar(0.5.*(box(2:end)+box(1:(end-1))),B)
plot(s,ys.*n.*dbox,'g')
scatter(estmgs,interp1(s,ys.*n.*dbox,estmgs),[],'rd','filled')
plot([s(1),s(end)],[lvs95.*n.*dbox,lvs95.*n.*dbox],'b:')
plot([s(1),s(end)],[lvs68.*n.*dbox,lvs68.*n.*dbox],'b--')
title('Within-sex CV')
xlabel('within-sex CV [%]')
axis([ls2,us2,0,1.05.*max(B)])
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
vals=[RS(1,:);Rkut]';
[~,~,binw]=histcounts(vals(:,1),frmw);
[~,~,binr]=histcounts(vals(:,2),frmr);
frame=zeros(length(frmr)-1,length(frmw)-1);
for i=1:(length(frmr)-1)
    for j=1:(length(frmw)-1)
        frame(i,j)=sum(binw==j&binr==(length(frmr)-i));
    end
end
set(gca,'FontSize',ftsz)
imagesc(frame,'XData',frmw(1:end-1)+0.5.*ds,'YData',frmr(end-1:-1:1)+0.5.*dk)
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