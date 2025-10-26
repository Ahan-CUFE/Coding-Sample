clear

result_cv=[];
result_hbic=[];
for end_year = 1994:2016
    temp = recursive_test(end_year);
    result_dl0 = temp.dl0;
    result_hbic1 = temp.hbic;
    result_cv  = [result_cv  result_dl0];
    result_hbic = [result_hbic result_hbic1];
end

function result  =  recursive_test(end_year)
%% all data loaded
allfactors = csvread('factors.csv',1,0);
date = allfactors(:,1);
rf = allfactors(:,2);
factors = allfactors(:,3:end);

L = length(date);
P = size(factors,2);

% test portfolios

port_3x2 = csvread('port_3x2.csv',0,1);
port_3x2 = port_3x2 - rf*ones(1,size(port_3x2,2));

% other information
summary = readtable('summary.csv');
factorname = summary.Row;
factorname_full = summary.Descpription;
year_pub = summary.Year;
year_end = summary.Year_end;

port_3x2_id = readtable('port_3x2_id.csv');


%% form a smaller set of portfolios for bivariate sorted porfolios
kk = 10; % minimun number of stocks in a portfolio
% 此处选出最小股票数大于10的以及截止时间之前的。

include_3x2 = find(port_3x2_id.min_stk6>=kk)';
port_index = find(port_3x2_id.Year<end_year)';

port_3x2b = [];

for i = 1:P
    if ismember(i,include_3x2) && ismember(i,port_index)
        port_3x2b = [port_3x2b port_3x2(:,(i*6-5):(i*6))];
    end
end

%%以下将在port3x2b中进行递归测试

%找出在end_year 之前的因子 summary.Year
factor_index = find(summary.Year<end_year);

% for i  = 1:length(port_index)
% 从1977年6月开始， 将时间维度取出
port_3x2t = port_3x2b(1:6+12*(end_year-1977),:);

Ri = port_3x2t';
ControlFactor = factors(:,factor_index);
Testlist = find(summary.Year == end_year);
TestFactor = factors(:,Testlist);
%%

% [size(port_3x2b,2) , size(ControlFactor,2) ,size(TestFactor,2)]

gt = TestFactor(1:6+12*(end_year-1977),:)'; % test factor
ht = ControlFactor(1:6+12*(end_year-1977),:)'; % control factor

result_dl0=[];
result_hbic=[];

for j = 1:size(gt,1)
    disp(j)
    seed_num = 200;
    gt1 = gt(j,:); % test factor

    %% CV
    N = size(Ri,1);
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
    model_Dl0_cv = DESDAR(Ri,gt1,ht,opts1,Kfld,seed_num,sim_times);
    model_Dl0.Year = end_year;
    model_Dl0.factorid = Testlist(j);
    model_Dl0.factorname = factorname(Testlist(j));
    model_Dl0.factornamefull = factorname_full(Testlist(j));
    model_Dl0.controlnum = size(ht,1);
    model_Dl0.assetnum = size(Ri,1);
    model_Dl0.tstat_Dl0= model_Dl0_cv.tstat_Dl0_final;
    model_Dl0.gamma = model_Dl0_cv.gamma_final;
    result_dl0 = [result_dl0 model_Dl0];


    %% HBIC

    N = size(Ri,1);
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
    model_Dl0_hbic = DESDAR(Ri,gt1,ht,opts2,Kfld,seed_num,sim_times);

    model_Dl01.Year = end_year;
    model_Dl01.factorid = Testlist(j);
    model_Dl01.factorname = factorname(Testlist(j));
    model_Dl01.factornamefull = factorname_full(Testlist(j));
    model_Dl01.controlnum = size(ht,1);
    model_Dl01.assetnum = size(Ri,1);
    model_Dl01.tstat_Dl0_hbic = model_Dl0_hbic.tstat_Dl0_final;
    model_Dl01.gamma = model_Dl0_hbic.gamma_final;
    result_hbic = [result_hbic model_Dl01];
    
end

% result = [result_hbic result_dl0];
    result.dl0 = result_dl0;
    result.hbic = result_hbic;

end

