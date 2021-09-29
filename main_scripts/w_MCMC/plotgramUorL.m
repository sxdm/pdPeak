function[]=plotgramUorL(data,sex,frcsd,f1)
if any(isnan(sex))
    sext=data<mean(data);
else
    sext=sex;
end
        phy=erf(frcsd/sqrt(2));
        avsd=mean([std(data(sext==0)),std(data(sext==1))]);
        trin=phy.*mean([sum(sext==0),sum(sext==1)]);
        l=2.*frcsd.*avsd;
        trin2=sqrt(2.*trin+0.25)-0.5;
        r=l./(2.*(trin2-1));
        width=r.*0.5;%
n=length(data);
tn=10000;
st=inf;
rdi=(1:n)';
for i=1:tn
    rdit=randperm(n)';
    y=plotgram(data(rdit),width);
    sy=sum(y);
    if sy<st
        st=sy;
        rdi=rdit;
    end
end
y=plotgram(data(rdi),width);
data=data(rdi);
sex=sex(rdi);
figure(f1)
clf
title('Data plot')
hold on
for i=1:length(data)
    if sex(i)==0
        rectangle('Position',[data(i)-width,y(i)-width,2.*width,2.*width],...
            'Curvature',[1,1],'LineStyle','non','FaceColor','b')
    elseif sex(i)==1
        rectangle('Position',[data(i)-width,y(i)-width,2.*width,2.*width],...
            'Curvature',[1,1],'LineStyle','non','FaceColor','r')
    else
        rectangle('Position',[data(i)-width,y(i)-width,2.*width,2.*width],...
            'Curvature',[1,1],'LineStyle','non','FaceColor',[0.3,0.3,0.3])
    end
end
hold off
axis equal
axis([min(data)-2.*width,max(data)+2.*width,0,max(y)+2.*width])
end
function[y]=plotgram(data,width)
y=repmat(-width,size(data));
yc=y;
for i=1:length(data)
    for j=1:length(data)
        if abs(data(j)-data(i))>=2*width
            yc(j)=width;
        else
            yc(j)=y(j)+sqrt(4.*width.^2-(data(j)-data(i)).^2);
        end
    end
       y(i)=max(yc);
end
end