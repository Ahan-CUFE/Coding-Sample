function result  = Postlaso_hbic_forsimu(Ri,gt,ht,opts,Kfld,seednum,sim_times)

if isempty(seednum)
    seednum = 100;
end
rng(seednum);
% dim of gt
d = size(gt,1);
sel2 = [];
sel3 = [];

% 1st and 2nd selection
for i = 1:d
    TSCVout =  D_TSCV(Ri, gt(i,:), ht,Kfld,opts,1,seednum);
    sel1 = TSCVout.sel1;
    sel2 = [sel2; TSCVout.sel2];
    sel3 = [sel3; TSCVout.sel3];
end
sel1 = unique(sel1);
sel2 = unique(sel2);
sel3 = unique(sel3);

% post-selection estimation and inference
dsout = infer(Ri, gt, ht, sel1, sel2, sel3);
ssout = infer(Ri, gt, ht, sel1, [], sel3);

% output for Double Selection
result.lambdag_ds = dsout.lambdag;
result.se_ds = dsout.se;
result.gamma_ds = dsout.gamma;

% output for Single Selection
result.lambdag_ss = ssout.lambdag;
result.se_ss = ssout.se;
result.gamma_ss = ssout.gamma;

% % selection results
result.sel1 = sel1;
result.sel2 = sel2;
result.sel3 = sel3;
select1 = unique([sel1;sel2]);
result.select = select1;

end


%%
% The function for cross-validation over time
%

function output = D_TSCV(Ri, gt, ht, Kfld, opts , Jrep,seednum)

if isempty(seednum)
    seednum = 100;
end

% data information
[p,T] = size(ht);
n = size(Ri,1);
alpha = opts.alpha;
lambda = opts.lambda0;

tmp3 = nancov([ht; Ri]');
Ch0   = tmp3((p+1):end,1:p);
tmp3b = nancov([gt; Ri]');
Cg0 = tmp3b(2:end,1);
ER0 = nanmean(Ri,2);

beta = NaN(n,p);
for i = 1:p
    beta(:,i) = Ch0(:,i)/nanvar(ht(i,:));
end
penalty = mean(beta.^2,1);
penalty = penalty./mean(penalty); % normalize the level

nomissing = (sum(isnan([ht;gt]),1)==0)';
%%

opts1 = struct('standardize',false,'lambda',lambda,'alpha',alpha);
model1 = glmnet(Ch0*(diag(penalty)), ER0,'gaussian',opts1);
coff1 = model1.beta;

[ebeta1, lam1, ~, ~, ~] = postLassoHBIC(coff1, Ch0*(diag(penalty)),  ER0, lambda, true);

opts2 = struct('standardize',false,'lambda',lambda,'alpha',alpha);
model2 = glmnet(Ch0*(diag(penalty)), Cg0,'gaussian',opts2);
coff2 = model2.beta;
[ebeta2, lam2, ~, ~, ~] = postLassoHBIC(coff2, Ch0*diag(penalty), Cg0, lambda, true);


opts3 = struct('intr',false,'standardize',true,'lambda',lambda,'alpha',alpha);
model3 = glmnet(ht(:,nomissing)', gt(:,nomissing)','gaussian',opts3);
coff3 = model3.beta;

[ebeta3, lam3, ~, ~, ~] = postLassoHBIC(coff3, ht(:,nomissing)', gt(:,nomissing)', lambda, false);




% 

sel1 = find(ebeta1 ~= 0);
sel2 = find(ebeta2 ~= 0);
sel3 = find(ebeta3 ~= 0);
output.sel1 = sel1;
output.sel2 = sel2;
output.sel3 = sel3;
output.lam1 = lam1;
output.lam2 = lam2;
output.lam3 = lam3;
end


%%
% the function for estimation and inference
%

function output = infer(Ri, gt, ht, sel1, sel2, sel3)

n = size(Ri,1);
p = size(ht,1);
d = size(gt,1);

tmp1  = nancov([gt',Ri']);
cov_g = tmp1((d+1):end,1:d);
tmp2  = nancov([ht',Ri']);
cov_h = tmp2((p+1):end,1:p);

ER    = mean(Ri,2);
%
M0 = eye(n) - ones(n,1)*inv(ones(n,1)'*ones(n,1))*ones(n,1)';

nomissing = find(sum(isnan([ht;gt]),1)==0);
Lnm = length(nomissing);
select = unique([sel1;sel2]);


X = [cov_g, cov_h(:,select)];
lambda_full = inv(X'*M0*X)*(X'*M0*ER);
lambdag = lambda_full(1:d);
clear X

% For double selection inference: AVAR
zthat = NaN(d,Lnm);
for i = 1:d
    M_mdl = eye(Lnm) - ht(sel3,nomissing)'*inv(ht(sel3,nomissing)*ht(sel3,nomissing)')*ht(sel3,nomissing);
    zthat(i,:) = M_mdl*gt(i,nomissing)';
    clear M_mdl
end
Sigmazhat = zthat*zthat'/Lnm;

temp2  =  0;
ii = 0;
for l = nomissing
    ii = ii+1;
    mt = 1-lambda_full'*[gt(:,l);ht(select,l)];
    temp2 = temp2 + mt^2*(inv(Sigmazhat)*zthat(:,ii)*zthat(:,ii)'*inv(Sigmazhat));
end

avar_lambdag = diag(temp2)/Lnm;
se = sqrt(avar_lambdag/Lnm);
clear temp2

% scaled lambda for DS
vt = [gt(:,nomissing);ht(select,nomissing)];
V_bar = vt - mean(vt,2)*ones(1,Lnm);
var_v = V_bar*V_bar'/Lnm;
gamma = diag(var_v).*lambda_full;
clear X vt V_bar var_v lambda_full

output.lambdag = lambdag;
output.se = se;
output.gamma = gamma;

end



function [betaBest, lambdaBest, betaPost, intercept, hbicVal] = postLassoHBIC(coff, X, Y, lambdaSeq, includeIntercept)
% postLassoHBIC   Perform Post-Lasso OLS with HBIC selection
% 
%   [betaBest, lambdaBest, betaPost, intercept, hbicVal] = postLassoHBIC(
%       coff, X, Y, lambdaSeq, includeIntercept)
% 
% Inputs:
%   coff            - p x L matrix of Lasso coefficients
%   X               - n x p design matrix
%   Y               - n x 1 response vector
%   lambdaSeq       - 1 x L vector of lambda values
%   includeIntercept- logical; if true, fit intercept in OLS
% 
% Outputs:
%   betaBest        - p x 1 Post-Lasso coefficients at best lambda
%   lambdaBest      - scalar; best lambda value (HBIC-minimizing)
%   betaPost        - p x L matrix of Post-Lasso coefficients for each lambda
%   intercept       - L x 1 vector of intercepts (or empty if includeIntercept=false)
%   hbicVal         - 1 x L vector of HBIC values

[n, p] = size(X);
L = size(coff,2);

betaPost = zeros(p, L);
if includeIntercept
    intercept = zeros(L,1);
else
    intercept = [];
end

hbicVal = nan(1,L);

for j = 1:L
    sup = find(coff(:,j) ~= 0);
    k   = numel(sup);
    
    if k > 0
        if includeIntercept
            Xsub = [ones(n,1), X(:,sup)];
            b    = Xsub \ Y;
            intercept(j)      = b(1);
            betaPost(sup,j)   = b(2:end);
            yhat               = Xsub * b;
        else
            Xsub = X(:,sup);
            b    = Xsub \ Y;
            betaPost(sup,j)   = b;
            yhat               = Xsub * b;
        end
    else
        if includeIntercept
            intercept(j) = mean(Y);
            yhat         = intercept(j)*ones(n,1);
        else
            yhat         = zeros(n,1);
        end
    end
    
    residual = Y - yhat;
    sigmaEst = sqrt(mean(residual.^2));
    hbicVal(j) = log(sigmaEst) + (log(log(n)) * log(p) / n) * k;
end

% pick best
[~, idx]   = min(hbicVal);
lambdaBest = lambdaSeq(idx);
betaBest   = betaPost(:, idx);

end