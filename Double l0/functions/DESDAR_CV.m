function result  = DESDAR_CV(Ri,gt,ht,opts,Kfld,seednum)
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
L = length(opts.tau);

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

cvm1 = NaN(L,Kfld,Jrep);
cvm2 = NaN(L,Kfld,Jrep);
cvm3 = NaN(L,Kfld,Jrep);

cvm11 = [];
cvm22 = [];
cvm33 = [];


nomissing = (sum(isnan([ht;gt]),1)==0)';
if opts.isCV == 1   
    for j = 1:Jrep
    
        rng(seednum+j)
        indices = crossvalind('Kfold',T,Kfld);
    
        for k = 1:Kfld
    
            % divide the train and test samples
            test = (indices == k);
            train = (indices ~= k);
    
            Ri_train = Ri(:,train);
            ht_train = ht(:,train);
            gt_train = gt(:,train);
            tmp1 = nancov([ht_train; Ri_train]');
            Ch_train = tmp1((p+1):end,1:p);
            tmp1b = nancov([gt_train; Ri_train]');
            Cg_train = tmp1b(2:end,1);
            ER_train = nanmean(Ri_train,2);
    
            ht_train = ht(:,train & nomissing);
            gt_train = gt(:,train & nomissing);
    
            Ri_test = Ri(:,test);
            ht_test = ht(:,test);
            gt_test = gt(:,test);
            tmp2 = nancov([ht_test; Ri_test]');
            Ch_test = tmp2((p+1):end,1:p);
            tmp2b = nancov([gt_test; Ri_test]');
            Cg_test = tmp2b(2:end,1);
            ER_test = nanmean(Ri_test,2);
    
 
    
            
            coff1=[];
            coff2=[];
            coff3=[];


            opts1.tau0 = opts.tau0;
            opts1.mu1  = opts.mu1;
            opts1.outputmode = opts.outputmode;
            opts1.J = opts.J;
            opts1.tau =opts.tau;
            opts1.delta = opts.delta;
            opts1.ifintercept = opts.ifintercept;
            opts1.initial = opts.initial;
            model1 = AESDAR(Ch_train*(diag(penalty)),ER_train,opts1);
            coff1 = model1.solutionpath;
            intercept_1 = model1.intercept;


            opts2.tau0 = opts.tau0;
            opts2.mu1  =  opts.mu1;
            opts2.outputmode = opts.outputmode;
            opts2.J = opts.J;
            opts2.tau =opts.tau;
            opts2.delta = opts.delta;
            opts2.ifintercept = opts.ifintercept;
            opts2.initial = opts.initial;
            model2 = AESDAR(Ch_train*diag(penalty),Cg_train,opts2);
            coff2 = model2.solutionpath;
            intercept_2= model2.intercept;


            opts3.tau0 = opts.tau0;
            N3 = size(ht_train,2);
            opts3.mu1 = (size(ht_train,1))/N3;
            opts3.outputmode = opts.outputmode;
            opts3.J = opts.J;
            opts3.tau =opts.tau;
            opts3.delta = opts.delta;
            opts3.ifintercept=0;
            opts3.initial = opts.initial;
            model3 = AESDAR((ht_train)',(gt_train)',opts3);
            coff3 = model3.solutionpath;

            ER_pred = Ch_test*diag(penalty)*coff1+intercept_1;
            Cg_pred = Ch_test*diag(penalty)*coff2+intercept_2;
            gt_pred = ht_test'*coff3;


            LL1 = size(coff1,2);
            LL2 = size(coff2,2);
            LL3 = size(coff3,2);

            cvm1(1:LL1,k,j) = mean((repmat(ER_test,1,LL1) - ER_pred).^2,1)';
            cvm2(1:LL2,k,j) = mean((repmat(Cg_test,1,LL2) - Cg_pred).^2,1)';
            cvm3(1:LL3,k,j) = nanmean((repmat(gt_test',1,LL3) - gt_pred).^2,1)';
        end
    
        cvm11 = [cvm11, cvm1(:,:,j)];
        cvm22 = [cvm22, cvm2(:,:,j)];
        cvm33 = [cvm33, cvm3(:,:,j)];
    end
    
    cv_sd1 = std(cvm11',1)/sqrt(Kfld*Jrep);
    cv_sd2 = std(cvm22',1)/sqrt(Kfld*Jrep);
    cv_sd3 = std(cvm33',1)/sqrt(Kfld*Jrep);
    
    cvm111 = mean(cvm11,2);
    cvm222 = mean(cvm22,2);
    cvm333 = mean(cvm33,2);
    
    [~,l_sel1] = min(cvm111);
    [~,l_sel2] = min(cvm222);
    [~,l_sel3] = min(cvm333);
    


%   normalize the X 
    [sX,D]  =  normalize((Ch0-mean(Ch0))*(diag(penalty)));
    optse1.T = l_sel1*opts.tau0;                 % number of features to be extract each step 
    optse1.alpha = 0;                  % a parameter for numerical stable 
    optse1.J = opts.J;              % maxiteration number in SDCA  
    optse1.initial = zeros(p,1);
    optse1.tau =opts.tau;
    optse1.ifintercept = opts.ifintercept;
    % run SDAR
    [ebeta1, ed1, Ac1, nIter1] = esdar(sX,ER0-mean(ER0),optse1);

    optse2.T = l_sel2*opts.tau0;                 % number of features to be extract each step 
    optse2.alpha = 0;                  % a parameter for numerical stable 
    optse2.J = opts.J;              % maxiteration number in SDCA  
    optse2.initial = zeros(p,1);
    optse2.tau =opts.tau;
    optse2.ifintercept=opts.ifintercept;
    % run SDAR
    [ebeta2, ed2, Ac2, nIter2] = esdar(sX,Cg0-mean(Cg0),optse2);
    
    [sX2,D]  =  normalize((ht(:,nomissing))');
    optse3.T = l_sel3*opts.tau0;                 % number of features to be extract each step 
    optse3.alpha = 0;                  % a parameter for numerical stable 
    optse3.J = opts.J;              % maxiteration number in SDCA  
    optse3.initial = zeros(p,1);
    optse3.tau =opts.tau;
    optse3.ifintercept =0;
    [ebeta3, ed3, Ac3, nIter3] = esdar(sX2,(gt(:,nomissing))',optse3);


else
  
    opts1.tau0 = opts.tau0;
    opts1.mu1  = opts.mu1;
    opts1.outputmode = opts.outputmode;
    opts1.J = opts.J;
    opts1.tau =opts.tau;
    opts1.delta = opts.delta;
    opts1.ifintercept =opts.ifintercept;
    opts1.initial = opts.initial;
    model1 = AESDAR(Ch0*(diag(penalty)),ER0,opts1);
    ebeta1 = model1.beta_hbic;

    opts2.tau0 = opts.tau0;
    opts2.mu1  =  opts.mu1;
    opts2.outputmode = opts.outputmode;
    opts2.J = opts.J;
    opts2.tau =opts.tau;
    opts2.delta = opts.delta;
    opts2.ifintercept = opts.ifintercept;
    opts2.initial = opts.initial;
    model2 = AESDAR(Ch0*(diag(penalty)),Cg0,opts2);
    ebeta2 = model2.beta_hbic;

    opts3.tau0 = opts.tau0;
    N3 = size(ht,2);
    opts3.mu1 = (size(ht,1))/N3;
    opts3.outputmode = opts.outputmode;
    opts3.J = opts.J;
    opts3.tau =opts.tau;
    opts3.delta = opts.delta;
    opts3.ifintercept = 0;
    opts3.initial = opts.initial;
    model3 = AESDAR((ht(:,nomissing))',(gt(:,nomissing))',opts3);
    ebeta3 = model3.beta_hbic;








end

sel1 = find(ebeta1 ~= 0);
sel2 = find(ebeta2 ~= 0);
sel3 = find(ebeta3 ~= 0);
output.sel1 = sel1;
output.sel2 = sel2;
output.sel3 = sel3;


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



