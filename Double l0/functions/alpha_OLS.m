
% Input:
% gt, esting factors
% ht, factor model

% output:
% output for the No Selection OLS Estimation
% alpha_ols, se_ols

function result = alpha_OLS(gt, control)

% data information


p = size(control,1);
nomissing = find(sum(isnan([control;gt]),1)==0);
n = length(nomissing);
%% For no selection OLS
% no selection estimate

% X_zero = [ones(n,1),control(:,nomissing)'];
% lambda_full_zero = inv(X_zero'*X_zero)*(X_zero'*gt(:,nomissing)');
% alpha = lambda_full_zero(1);
% 
% residuals = gt(:,nomissing)'-X_zero*lambda_full_zero;
% SE = sqrt(diag((1/(n-p)).*residuals'*residuals*inv(X_zero'*X_zero)));
% se_alpha = SE;

[~,se,coeff] = hac(control(:,nomissing)',gt(:,nomissing)');

alpha1 = coeff(1);
se_alpha1  =se(1);
t1  =alpha1/se_alpha1;


result.t_alpha = t1;
result.alpha = alpha1;






