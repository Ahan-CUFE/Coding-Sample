
% Main results

clear; clc;
close all;

addpath('/glmnet_matlab')
addpath('/functions')
addpath('../data')


% factor
allfactors = csvread('factors.csv',1,0);
date = allfactors(:,1);
rf = allfactors(:,2); % risk free rate
factors = allfactors(:,3:end);

L = length(date);
P = size(factors,2);


% other information
summary = readtable('summary.csv');
factorname = summary.Row;
factorname_full = summary.Descpription;
year_pub = summary.Year;
year_end = summary.Year_end;

port_3x2_id = readtable('port_3x2_id.csv');

profitability = {'sgr', 'cto', 'os', 'zs', 'ps', 'chempia', 'ms', 'rna', 'pm', 'ato', 'chatoia', 'chpmia', 'roaq', 'gma', 'RMW', 'HXZ_ROE'};
profitability_ind = find(ismember(factorname,profitability));

%% FF
mkt_ind = find(ismember(factorname,'MktRf'));
smb_ind = find(ismember(factorname,'SMB'));
hml_ind = find(ismember(factorname,'HML'));
umd_ind = find(ismember(factorname,'UMD'));
rmw_ind = find(ismember(factorname,'RMW'));
cma_ind = find(ismember(factorname,'CMA'));

%% SDAR HBIC
cfp_ia_ind = find(ismember(factorname,'cfp_ia'));
pchcurrat_ind = find(ismember(factorname,'pchcurrat'));
dsti_ind = find(ismember(factorname,'dsti'));
sp_ind = find(ismember(factorname,'sp'));
umd_ind = find(ismember(factorname,'UMD'));
noa_ind = find(ismember(factorname,'noa'));

%% SDAR     

cfp_ia_ind = find(ismember(factorname,'cfp_ia'));
mve_ia_ind = find(ismember(factorname,'mve_ia'));
umd_ind = find(ismember(factorname,'UMD'));
ill_ind = find(ismember(factorname,'ill'));
HXZ_ROE_ind = find(ismember(factorname,'HXZ_ROE'));
LTR_ind  = find(ismember(factorname,'LTR'));

%% post lasso

HXZ_ROE_ind = find(ismember(factorname,'HXZ_ROE'));
chatoia_ind = find(ismember(factorname,'chatoia'));
mve_ia_ind = find(ismember(factorname,'mve_ia'));
acc_ind = find(ismember(factorname,'acc'));
LIQ_PS_ind = find(ismember(factorname,'LIQ_PS'));
ill_ind = find(ismember(factorname,'ill'));

%% form a smaller set of portfolios for bivariate sorted porfolios
FF3 = factors(:,[mkt_ind,smb_ind,hml_ind])';
FF5 = factors(:,[mkt_ind,smb_ind,hml_ind,rmw_ind,cma_ind])';
FF6 = factors(:,[mkt_ind,smb_ind,hml_ind,rmw_ind,cma_ind,umd_ind])';
%%
Dl03_hbic_ind = [pchcurrat_ind,cfp_ia_ind,dsti_ind];
Dl05_hbic_ind = [pchcurrat_ind,cfp_ia_ind,dsti_ind,sp_ind,umd_ind];
Dl06_hbic_ind = [pchcurrat_ind,cfp_ia_ind,dsti_ind,sp_ind,umd_ind,noa_ind];
Dl03_hbic = factors(:,Dl03_hbic_ind)';
Dl05_hbic = factors(:,Dl05_hbic_ind)';
Dl06_hbic = factors(:,Dl06_hbic_ind)';
%%
Dl06_ind = [mve_ia_ind,cfp_ia_ind,umd_ind,ill_ind,HXZ_ROE_ind,LTR_ind];
Dl06 = factors(:,Dl06_ind)';
%%
DPL_ind6 = [HXZ_ROE_ind,chatoia_ind,mve_ia_ind,acc_ind,LIQ_PS_ind,ill_ind];
DPL6 = factors(:,DPL_ind6)';
%%



result=[];
for i = 1:length(profitability)
    disp(i)
    gt = factors(:,profitability_ind(i))'; % test factor
    j=profitability_ind(i);

    if ismember(j, [mkt_ind, smb_ind, hml_ind]) % If the factor is in FF3 model
        alpha_FF3 = NaN;
        tstat_FF3 = NaN;
    else
        result_FF3 = alpha_OLS(gt, FF3);
        alpha_FF3  = result_FF3.alpha;
        tstat_FF3  = result_FF3.t_alpha;
    end

    if ismember(j, [mkt_ind, smb_ind, hml_ind, rmw_ind, cma_ind]) % If the factor is in FF5 model
        alpha_FF5 = NaN;
        tstat_FF5 = NaN;
    else
        result_FF5= alpha_OLS(gt, FF5);
        alpha_FF5  = result_FF5.alpha;
        tstat_FF5  = result_FF5.t_alpha;
    end
    
    if ismember(j, [mkt_ind, smb_ind, hml_ind, rmw_ind, cma_ind, umd_ind]) % If the factor is in FF6 model
        alpha_FF6 = NaN;
        tstat_FF6 = NaN;
    else
        result_FF6 = alpha_OLS(gt, FF6);
        alpha_FF6  = result_FF6.alpha;
        tstat_FF6  = result_FF6.t_alpha;
    end

%%  Double SDAR + HBIC
    if ismember(j, Dl03_hbic_ind) % If the factor is in Dl03 model
        alpha_Dl03_hbic = NaN;
        tstat_Dl03_hbic = NaN;
    else
    
        result_Dl03_hbic = alpha_OLS(gt, Dl03_hbic);
        alpha_Dl03_hbic  = result_Dl03_hbic.alpha;
        tstat_Dl03_hbic  = result_Dl03_hbic.t_alpha;
    end


    if ismember(j, Dl05_hbic_ind) % If the factor is in Dl05 model
        alpha_Dl05_hbic = NaN;
        tstat_Dl05_hbic = NaN;
    else

        result_Dl05_hbic = alpha_OLS(gt, Dl05_hbic);
        alpha_Dl05_hbic  = result_Dl05_hbic.alpha;
        tstat_Dl05_hbic  = result_Dl05_hbic.t_alpha;
    end

    if ismember(j, Dl06_hbic_ind) % If the factor is in Dl06 model
        alpha_Dl06_hbic = NaN;
        tstat_Dl06_hbic = NaN;
    else
        result_Dl06_hbic = alpha_OLS(gt, Dl06_hbic);
        alpha_Dl06_hbic  = result_Dl06_hbic.alpha;
        tstat_Dl06_hbic  = result_Dl06_hbic.t_alpha;
    end

%% Double SDAR
    if ismember(j, Dl06_ind) % If the factor is in Dl06 model
        alpha_Dl06 = NaN;
        tstat_Dl06 = NaN;
    else
        result_Dl06 = alpha_OLS(gt, Dl06);
        alpha_Dl06  = result_Dl06.alpha;
        tstat_Dl06  = result_Dl06.t_alpha;
    end

%% Double Post Lasso HBIC

    if ismember(j, DPL_ind6) % If the factor is in Dl06 model
        alpha_DPL6 = NaN;
        tstat_DPL6 = NaN;
    else
        result_Dl06 = alpha_OLS(gt, DPL6);
        alpha_DPL6  = result_Dl06.alpha;
        tstat_DPL6  = result_Dl06.t_alpha;
    end
%%




    temp = table(alpha_FF3,tstat_FF3,alpha_FF5,tstat_FF5,alpha_FF6,tstat_FF6, ...
      alpha_Dl03_hbic, tstat_Dl03_hbic,alpha_Dl05_hbic,tstat_Dl05_hbic,alpha_Dl06_hbic ,tstat_Dl06_hbic ...
      ,alpha_Dl06, tstat_Dl06, alpha_DPL6,tstat_DPL6 );

    result = [result ;temp];

end

result1 = abs(result);
result = [result ; nanmean(result1)];
result = [result ; sum(result1>1.96)];
result = [result ; sum(result1>2.8)];
writetable(result,'profitability.csv')