
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
umd_ind = find(ismember(factorname,'UMD'));
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

FF3 = factors(:,[mkt_ind,smb_ind,hml_ind])';
FF5 = factors(:,[mkt_ind,smb_ind,hml_ind,rmw_ind,cma_ind])';
FF6 = factors(:,[mkt_ind,smb_ind,hml_ind,rmw_ind,cma_ind,umd_ind])';




result=[];
result_hbic = [];
for j = 147:length(factorname_full)
    disp(j)
    gt = factors(:,j)'; % test factor
    ht = factors(:,setdiff(1:P, j))';  % control factor

%%
    N = size(Ri,2);
    opts1.tau=1;
    opts1.J = 400;
    opts1.tau0 = 3;
    opts1.isCV = 1;
    opts1.outputmode = 1;
    opts1.mu1 = (size(ht,1))/N;
    opts1.ifintercept = 1;
    opts1.initial='default';
    opts1.delta = 1e-6;
    Kfld = 10;
    sim_times = 200;
    model_Dl0 = DESDAR(Ri',gt,ht,opts1,Kfld,seed_num,sim_times);
    result = [result ; model_Dl0 ];




%%

    N = size(Ri,2);
    opts2.tau=1;
    opts2.J = 400;
    opts2.tau0 = 1;
    opts2.isCV = 0;
    opts2.outputmode = 2;
    opts2.mu1 =(size(ht,1))/N;
    opts2.delta = 1e-6;
    opts2.ifintercept = 1;
    opts2.initial=[];
    Kfld = 10;
    sim_times = 200;
    model_Dl0_hbic = DESDAR(Ri',gt,ht,opts2,Kfld,seed_num,sim_times);
    result_hbic = [ result_hbic ; model_Dl0_hbic];

    %%

    opts3.lambda0 = exp(linspace(0,-35,100));
    opts3.alpha = 1;
    Kfld = 10;
    sim_times = 200;
    model_DBS_hbic = Postlaso_hbic(Ri',gt,ht,opts3,Kfld,seed_num,sim_times);
    lambda_DBL_hbic = model_DBS_hbic.gamma_final;
    tstat_DBL_hbic = model_DBS_hbic.tstat_Dl0_final;
    lambda_sl_hbic = model_DBS_hbic.gamma_ss_final;
    tstat_sl_hbic = model_DBS_hbic.tstat_sl0_final;
    result = [result;model_DBS_hbic];




    %%

    % model_ols = PriceRisk_OLS(Ri', gt, ht);
    % tstat_ols = model_ols.lambdag_ols/model_ols.se_ols;
    % lambda_ols = model_ols.lambda_ols(1);
    % 
    % % only control FF3 by OLS
    % model_FF3 = PriceRisk_OLS(Ri', gt, FF3);
    % tstat_FF3 = model_FF3.lambdag_ols/model_FF3.se_ols;
    % lambda_FF3 = model_FF3.lambda_ols(1);
    % 
    % 
    % model_FF5 = PriceRisk_OLS(Ri', gt, FF5);
    % tstat_FF5 = model_FF5.lambdag_ols/model_FF5.se_ols;
    % lambda_FF5 = model_FF5.lambda_ols(1);
    % 
    % 
    % model_FF6 = PriceRisk_OLS(Ri', gt, FF6);
    % tstat_FF6 = model_FF6.lambdag_ols/model_FF6.se_ols;
    % lambda_FF6 = model_FF6.lambda_ols(1);
    % 
    % 
    % 
    % 
    % 
    % temp = table(lambda_ols*1e4,tstat_ols,lambda_FF3*1e4,tstat_FF3,lambda_FF5*1e4,tstat_FF5,lambda_FF6*1e4,tstat_FF6);
    % temp = table(lambda_ols*1e4,tstat_ols,lambda_FF5*1e4,tstat_FF5);
    % result = [result ;temp];
    
end







