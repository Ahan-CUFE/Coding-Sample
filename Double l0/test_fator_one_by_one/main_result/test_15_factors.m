
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

% test factor individually
for j = 1:length(TestList)
    disp(j)

    gt = TestFactor(:,j)'; % test factor
    ht = ControlFactor'; % control factor

    
    N = size(Ri,2);
    opts1.tau=1;
    opts1.J = 400;
    opts1.tau0 = 3;
    opts1.isCV = 1;
    opts1.outputmode = 1;
    opts1.mu1 =(size(ContrlList,1))/N;
    opts1.ifintercept = 1;
    opts1.initial='default';
    opts1.delta = 1e-6;
    Kfld = 10;
    sim_times = 200;
    model_Dl0 = DESDAR(Ri',gt,ht,opts1,Kfld,seed_num,sim_times);
    lambda_Dl0 = model_Dl0.gamma_final;
    tstat_Dl0 = model_Dl0.tstat_Dl0_final;
    lambda_sl0 = model_Dl0.gamma_ss_final;
    tstat_sl0 = model_Dl0.tstat_sl0_final;
    % 

% lambda_Dl0=1;
%    tstat_Dl0=1; 
% lambda_sl0=1;
% tstat_sl0=1;

    N = size(Ri,2);
    opts2.tau = 1;
    opts2.J = 400;
    opts2.tau0 = 1;
    opts2.isCV = 0;  % 1 CV 0 HBIC 
    opts2.outputmode = 2;  % 1 for rmse 2 for hbic
    opts2.mu1 =(size(ContrlList,1))/N;
    opts2.delta = 1e-6;
    opts2.ifintercept = 1;
    opts2.initial=[];
    Kfld = 10;
    sim_times = 200;
    model_Dl0_hbic = DESDAR(Ri',gt,ht,opts2,Kfld,seed_num,sim_times);
    
    tstat_Dl0_hbic = model_Dl0_hbic.tstat_Dl0_final;
    lambda_Dl0_hbic = model_Dl0_hbic.gamma_final;



%%
    opts3.lambda0 = exp(linspace(0,-35,100));
    opts3.alpha = 1;
    Kfld = 10;
    sim_times = 200;
    model_DBS_hbic = Postlaso_hbic(Ri',gt,ht,opts3,Kfld,seed_num,sim_times);
    lambda_DBL_hbic = model_DBS_hbic.gamma_final;
    tstat_DBL_hbic = model_DBS_hbic.tstat_Dl0_final;

    %%
    model_ds = DS(Ri', gt, ht,-log(tune_center(j,1)),-log(tune_center(j,2)),1,seed_num);
    tstat_ds = model_ds.lambdag_ds/model_ds.se_ds;
    lambda_ds = model_ds.gamma_ds(1);

    %%
    % Single-Selection results, replace with a huge tune2
    model_ss = DS(Ri', gt, ht,-log(tune_center(j,1)),-log(1),1,seed_num);
    tstat_ss = model_ss.lambdag_ds/model_ss.se_ds;
    lambda_ss = model_ss.gamma_ds(1);
%%
    % controlling everything, no selection, OLS
    model_ols = PriceRisk_OLS(Ri', gt, ht);
    tstat_ols = model_ols.lambdag_ols/model_ols.se_ols;
    lambda_ols = model_ols.lambda_ols(1);

    % only control FF3 by OLS
    model_FF3 = PriceRisk_OLS(Ri', gt, FF3);
    tstat_FF3 = model_FF3.lambdag_ols/model_FF3.se_ols;
    lambda_FF3 = model_FF3.lambda_ols(1);

    % time series average
    avg = nanmean(gt);
    tstat_avg = avg/nanstd(gt)*sqrt(sum(~isnan(gt)));

    % combine the results in a table
    temp = table(tstat_Dl0,lambda_Dl0,tstat_Dl0_hbic,lambda_Dl0_hbic,tstat_ds, lambda_ds,tstat_DBL_hbic,lambda_DBL_hbic,tstat_ss,lambda_ss,avg,tstat_avg,...
        lambda_ols,tstat_ols,lambda_FF3,tstat_FF3);

    result = [result;temp];

end

% extract outputs
lambda_Dl0 = result.lambda_Dl0*10000;
tstat_Dl0  = result.tstat_Dl0;

lambda_Dl0_hbic=result.lambda_Dl0_hbic*10000;
tstat_Dl0_hbic = result.tstat_Dl0_hbic;

lambda_ds = result.lambda_ds*10000; % bp
tstat_ds = result.tstat_ds;

lambda_DBL_hbic = result.lambda_DBL_hbic*10000;
tstat_DBL_hbic  = result.tstat_DBL_hbic;

lambda_ss = result.lambda_ss*10000; % bp
tstat_ss = result.tstat_ss;


lambda_FF3 = result.lambda_FF3*10000; % bp
tstat_FF3 = result.tstat_FF3;

avg = result.avg*10000; % bp
tstat_avg = result.tstat_avg;

lambda_ols = result.lambda_ols*10000; % bp
tstat_ols = result.tstat_ols;

% factor names for those factors since 2012
factornames = factorname_full(TestList);

result = table(TestList,factornames,lambda_Dl0,tstat_Dl0,lambda_Dl0_hbic,tstat_Dl0_hbic,lambda_ds,tstat_ds,lambda_DBL_hbic,tstat_DBL_hbic,lambda_ss,tstat_ss,lambda_FF3,...
    tstat_FF3,lambda_ols,tstat_ols,avg,tstat_avg);

% display the table
disp(result)