
% Main results

clear; clc;
close all;

addpath('/glmnet_matlab')
addpath('/functions')
addpath('../data')

% fix the random seed for any additional CV
seed_num = 200;

% factor
allfactors = csvread('factors.csv',1,0);
date = allfactors(:,1);
rf = allfactors(:,2); % risk free rate
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
FF3 = factors(:,[mkt_ind,smb_ind,hml_ind])';

% test factors since 2012
TestList = find(year_pub >= 2012);
TestFactor = factors(:,TestList);

result = [];
result_hbic = [];
% test factor individually
for j = 1:length(TestList)
    disp(j)

    gt = TestFactor(:,j)'; % test factor
    ht = ControlFactor'; % control factor
    
    opts1.lambda0 = exp(linspace(0,-35,100));
    opts1.alpha = 1;
    Kfld = 10;
    sim_times = 200;
    model_DBS_hbic = Postlaso_hbic(Ri',gt,ht,opts1,Kfld,seed_num,sim_times);
    lambda_DBL_hbic = model_DBS_hbic.gamma_final;
    tstat_DBL_hbic = model_DBS_hbic.tstat_Dl0_final;
    lambda_sl_hbic = model_DBS_hbic.gamma_ss_final;
    tstat_sl_hbic = model_DBS_hbic.tstat_sl0_final;
    result = [result;model_DBS_hbic];


% lambda_Dl0=1;
%    tstat_Dl0=1; 
% lambda_sl0=1;
% tstat_sl0=1;

%     N = size(Ri,2);
%     opts2.tau = 1;
%     opts2.J = 400;
%     opts2.tau0 = 1;
%     opts2.isCV = 0;  % 1 CV 0 HBIC 
%     opts2.outputmode = 2;  % 1 for rmse 2 for hbic
%     opts2.mu1 =(size(ContrlList,1))/N;
%     opts2.delta = 1e-6;
%     opts2.ifintercept = 1;
%     opts2.initial=[];
%     Kfld = 10;
%     sim_times = 200;
%     model_Dl0_hbic = DESDAR(Ri',gt,ht,opts2,Kfld,seed_num,sim_times);
% 
%     tstat_Dl0_hbic = model_Dl0_hbic.tstat_Dl0_final;
%     lambda_Dl0_hbic = model_Dl0_hbic.gamma_final;
% %     tstat_sl0_hbic = model_Dl0_hbic.tstat_sl0_final;
% %     lambda_sl0_hbic = model_Dl0_hbic.gamma_ss_final;
%     result_hbic = [result_hbic ; model_Dl0_hbic];

end
