% Main results

clear; clc;
close all;

addpath('/glmnet_matlab')
addpath('/functions')
addpath('../data')

% fix the random seed for any additional CV
seed_num = 100;

% factor
allfactors = csvread('factors.csv',1,0);
date = allfactors(:,1);
rf = allfactors(:,2);
factors = allfactors(:,3:end);

L = length(date);
P = size(factors,2);

% test portfolios
port_3x2 = csvread('port_3x2.csv',0,1);
port_3x2 = port_3x2 - rf*ones(1,size(port_3x2,2)); % excess return

% other information
summary = readtable('summary.csv');
factorname = summary.Row;
factorname_full = summary.Descpription;
year_pub = summary.Year;
year_end = summary.Year_end;

port_3x2_id = readtable('port_3x2_id.csv');

mkt_ind = find(ismember(factorname,'MktRf'));
smb_ind = find(ismember(factorname,'SMB'));
hml_ind = find(ismember(factorname,'HML'));
rmw_ind = find(ismember(factorname,'RMW'));
cma_ind = find(ismember(factorname,'CMA'));
%% form a smaller set of portfolios for bivariate sorted porfolios
kk = 10; % minimun number of stocks in a portfolio

include_3x2 = find(port_3x2_id.min_stk6>=kk)';
port_3x2b = [];

for i = 1:P
    if ismember(i,include_3x2)
        port_3x2b = [port_3x2b port_3x2(:,(i*6-5):(i*6))];
    end
end

Ri = port_3x2b; % test asset

% load tune_center for 200 randome seeds selected by cross-validations
% For each randome seed in 1:200, we run a cross-validation and find the
% best tuning parameter.
% Then we take the average for the 200 selected tuning parameters at the
% log scale as tune center, because we plot heat maps on log(lambda).
load tune_main.mat

% choose control factors before 2012
ContrlList = find(year_pub < 2012);
ControlFactor = factors(:,ContrlList);
FF5 = factors(:,[mkt_ind,smb_ind,hml_ind,rmw_ind,cma_ind])';
% test factors since 2012
TestList = find(year_pub >= 2012);
TestFactor = factors(:,TestList);

Useless_index = setdiff([1:150], [mkt_ind,smb_ind,hml_ind,rmw_ind,cma_ind])
Uselessfactor = factors(:,Useless_index)'

[n,T] = size(Ri);
d = size(Uselessfactor,1);
p = size(FF5,1);


ht = FF5;
gt = Uselessfactor;
Ri = Ri';

cov_h = NaN(n,p);
for nn = 1:n
    temp = nancov([Ri(nn,:)' ht']);
    cov_h(nn,:) = temp(1,2:end);
end
cov_g = NaN(n,d);
for nn = 1:n
    temp = nancov([Ri(nn,:)' gt']);
    cov_g(nn,:) = temp(1,2:end);
end






