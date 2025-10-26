function result = AESDAR(X,y,opts_ae)
%% normlize the data
if opts_ae.ifintercept
    meanX = mean(X);
    meany = mean(y);

else
    meanX = zeros(1,size(X,2));
    meany = zeros(1,size(y,2));
end

cX = X-meanX;
cy = y-meany;

[sX,D]  =  normalize(cX);  % X is n x P;
[n,p]   =  size(cX);
%
% sX = X;
% [n,p]   =  size(X);
% D = eye(135);

noiselev = 0;  % if can not know exactly just set small enough
                          % then the baseline ASDCA + BIC will work well 
%% Set parameters in ASDCA 
delta = opts_ae.delta;      % stop when norm(\beta(tau*k)-beta(tau*(k+1)))<delta
                   % if not known how to set just set small enough 
                   % then the baseline ASDCA + BIC will work well 
tau0 = opts_ae.tau0;       % adaptive scale, tau small accurate but may slower 
mu1 = opts_ae.mu1;         % stop when tau*k > mu*n 控制最大迭代步数 

%% Run ASDAR 
k = 0; 
residual = 0;
change = delta  + 1;
% store the solutions at each tau*k in Solutionpth and plot solution path  
% and use BIC select one for safeguard
kmax = floor(mu1*n/tau0); %最大迭代步数
Solutionpath = zeros(p,kmax); 
interceptpath = zeros(1,kmax);
Sresidual = [];
Siter = [];
% output by outputmode,  1 for MSE; 2 for HBIC 
outputmode = opts_ae.outputmode;  
HBic = [];
betaold = zeros(p,1);


while k < kmax &&  change > delta
    k = k + 1;
    %% Set parameters in SDAR 

    opts.T = tau0*k;                 % number of features to be extract each step 
    opts.alpha = 0;                  % a parameter for numerical stable 
    opts.J = opts_ae.J;              % maxiteration number in SDCA  
    opts.tau =opts_ae.tau;

    if isempty(opts_ae.initial)      % warmstart for HBIC+SDAR 
        Xtrans = X';
        gVec   = sum(X.^2,1);
        bk = zeros(p,1);
        dk = (Xtrans*(y - X*bk - meany)) ./ gVec';
        Hval = sqrt(gVec') .* abs(0.5*bk + 0.5*dk); %阈值算子
        sortedH = sort(Hval, 'descend');
        threshold = sortedH(tau0*k);
        activeSet = find(Hval >= threshold);
        fitRidge  = ridge_fun(sX(:,activeSet), cy, 0, opts_ae.ifintercept);
        bk(activeSet) = fitRidge(2:end); %系数项
        opts.initial = bk;
    else
        opts.initial =zeros(p,1);
    end


    % run SDAR
    [ebeta, ed, Ac, nIter] = esdar(sX,cy,opts);
    % opts.initial = ebeta;    % warm start 
    change = norm(D*(ebeta - betaold));
    betaold = ebeta;
    
    intercept = meany - meanX*D*ebeta;
    residual  = y-X*D*ebeta-intercept;
    Sresidual = [Sresidual, residual];
    Siter = [Siter, nIter];
    Solutionpath(:,k) = D*ebeta;
    interceptpath(:,k) = intercept;
    sigmaEst  = sqrt(mean((residual).^2));
    hbicVal   = log(sigmaEst) + log(log(n))*log(p)/n*tau0*k;
    HBic = [HBic,hbicVal];  
end   

%% Output  a solution 
if outputmode == 1
    % end
    % id  = find(Sresidual == min(Sresidual));
    result.solutionpath = Solutionpath;
    result.intercept = interceptpath;
    disp('solution path achieved')
else
    
    result.solutionpath = Solutionpath;
    result.intercept = interceptpath;
    id = find(HBic == min(HBic));
    betahat = Solutionpath(:,id(1));
    result.beta_hbic = betahat;
    result.HBic = HBic;
    disp('BIC mode is adopted')
end

end

function esti = ridge_fun(Xr, Yr, alpha_r, intcp)
    [nr, pr] = size(Xr);
    if intcp
        mx = mean(Xr,1);
        my = mean(Yr);
    else
        mx = zeros(1, pr);
        my = 0;
    end
    % Xc_r = bsxfun(@minus, Xr, mx);
    Xc_r =  Xr - mx ;
    Yc_r = Yr - my;
    coefs = (Xc_r'*Xc_r + alpha_r*eye(pr)) \ (Xc_r'*Yc_r);
    esti  = [my - mx*coefs; coefs];  % (pr+1)x1
end


