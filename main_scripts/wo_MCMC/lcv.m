function[out]=lcv(data)
%Plavcan, 1994, Am J Phys Anthropol 94.
slope=0.0214;intercept=-0.047;
x=std(data)./mean(data).*100;
out=exp(intercept+slope.*x);
end