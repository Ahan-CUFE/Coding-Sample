
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
port_202 = csvread('port202.csv',0,1);
port_202 = port_202 - rf*ones(1,size(port_202,2)); % excess return

% other information
summary = readtable('summary.csv');
factorname = summary.Row;
factorname_full = summary.Descpription;
year_pub = summary.Year;
year_end = summary.Year_end;
% 

Ri = port_202; % test asset

% choose control factors before 2012
ContrlList = find(year_pub < 2012);
ControlFactor = factors(:,ContrlList);

% test factors since 2012
TestList = find(year_pub >= 2012);
TestFactor = factors(:,TestList);

result = [];
result_DPL = [];
result_hbic = [];
% test factor individually
for j = 1:length(TestList)
    disp(j)

    gt = TestFactor(:,j)'; % test factor
    ht = ControlFactor'; % control factor

    N = size(Ri,2);
    %%
    % opts1.tau=1;
    % opts1.J = 400;
    % opts1.tau0 = 3;
    % opts1.isCV = 1;
    % opts1.outputmode = 1;
    % opts1.mu1 =(size(ht,1))/N;
    % opts1.ifintercept = 1;
    % opts1.initial='default';
    % opts1.delta = 1e-6;
    % Kfld = 10;
    % sim_times = 200;
    % model_Dl0 = DESDAR(Ri',gt,ht,opts1,Kfld,seed_num,sim_times);
    % result = [result ; model_Dl0 ];
    % 
    % %%
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
    result_hbic = [result_hbic; model_Dl0_hbic];

%     opts3.lambda0 = exp(linspace(0,-35,100));
%     opts3.alpha = 1;
%     Kfld = 10;
%     sim_times = 200;
%     model_DBS_hbic = Postlaso_hbic(Ri',gt,ht,opts3,Kfld,seed_num,sim_times);
%     lambda_DBL_hbic = model_DBS_hbic.gamma_final;
%     tstat_DBL_hbic = model_DBS_hbic.tstat_Dl0_final;
%     lambda_sl_hbic = model_DBS_hbic.gamma_ss_final;
%     tstat_sl_hbic = model_DBS_hbic.tstat_sl0_final;
%     result_DPL = [result_DPL;model_DBS_hbic];

end








