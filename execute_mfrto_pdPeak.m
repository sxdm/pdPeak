clear
close all
Nt=10;%if N>=Nt program uses MCMC, if N<Nt it does not use MCMC.

load metricData.mat

clearvars -except data Nt
if ~isvector(data)
    error('data should be a vector')
end
data=data(data>0);
data=data(:);
curdir=cd;
if length(unique(data))<4
    error('the number of unique data has to be >=4')
elseif length(data)>=Nt
    save('./main_scripts/w_MCMC/tmp.mat')
    cd ./main_scripts/w_MCMC
    executeMCMC
    executeMAP
    statistics
    load tmp.mat
    cd(curdir)
elseif length(data)<Nt
    save('./main_scripts/wo_MCMC/tmp.mat')
    cd ./main_scripts/wo_MCMC
    executeMAP
    statistics
    load tmp.mat
    cd(curdir)
end