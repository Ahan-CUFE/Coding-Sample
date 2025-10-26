clear;



N = 400;
T = 240;
d = 3;
p1= 4;
p2= 56;



lam_OLS_1 = [];%检验数据,人工选择初始变量回归
lam_OLS_2 = [];%OLS
lam_DBL   = [];%double lasso 
lam_SL    = [];%single lasso 
lam_Dl0   = [];
sel_Dl0   = [];
lam_OLS_cov=[];
lam_SL0    =[];
lam_Dl0_hbic=[];
lam_SL0_hbic=[];
sel_Dl0_hbic=[];

t_OLS_1 = [];%检验数据,人工选择初始变量回归
t_OLS_2 = [];%OLS
t_DBL   = [];%double lasso 
t_SL    = [];%single lasso 
t_DL0   = [];
t_SL0   = [];
t_bias_DBL=[];
t_bias_Dl0 =[];
t_bias_sl0 = [];

t_DL0_hbic=[];
t_SL0_hbic=[];
t_bias_Dl0_hbic=[];
t_bias_sl0_hbic=[];




Sim_time = 1000;

[covsigmaz,sigmah,eta,sigma_ut,betag,betah,ER,lambdag,chi1,betag1,betah1  ]=gendata_sima(N,T,d,p1,p2);



reccoff_DBL =zeros(Sim_time,p1+p2);
reccoff_DBL1=zeros(Sim_time,p1+p2);
reccoff_SL =zeros(Sim_time,p1+p2);
reccoff_Dl0=zeros(Sim_time,p1+p2);
reccoff_Dl01=zeros(Sim_time,p1+p2);
reccoff_Dl02=zeros(Sim_time,p1+p2);
reccoff_Dl03=zeros(Sim_time,p1+p2);
reccoff_Dl0_hbic=zeros(Sim_time,p1+p2);
reccoff_Dl01_hbic=zeros(Sim_time,p1+p2);
reccoff_Dl02_hbic=zeros(Sim_time,p1+p2);
reccoff_Dl03_hbic=zeros(Sim_time,p1+p2);


for i = 1:Sim_time

    display(i);
    seednum = 200+i;
    rng(seednum);
    ztht = mvnrnd(zeros(p1+p2+d,1),blkdiag(covsigmaz, sigmah),T);
    zt = ztht(:,1:d)';
    ht= ztht(:,d+1:end)';
    ht1=ht(1:4,:);
    gt = eta*ht+zt;
    nu =5;
    utz = mvnrnd(zeros(1, N),sigma_ut*diag(ones(N,1)), T); % 多元正态分布
    W = chi2rnd(nu, T, 1); % 自由度为 nu 的卡方分布
    ut = sqrt(nu ./ W) .* utz; % 每一行是一个 t 分布样本
    Ri= ER +betag*gt + betah*ht+ut';
%     Ri =ER+betag1*gt(1,:) + betah1*ht1+ut';

    ER2 = mean(Ri,2);
    % gt1 = gt(1,:);
    gt1 = gt;
    ht1 = ht(1:4,:);
    tmp1 = cov([gt1;Ri]');
    cov_g = tmp1(d+1:end,1:d);
    tmp2 = cov([ht1;Ri]');
    cov_h = tmp2(5:end,1:4);
    tmp3 = cov([ht;Ri]');
    cov_h2= tmp3(p1+p2+1:end,1:p1+p2);

%%  检验数据,人工选择初始变量回归
    X = [ones(N,1) cov_g(:,1) cov_h];
    M_1  =  eye(N) -  X*inv(X'*X)*X';
    P_1 = inv(X'*X)*X';
    lam1 = P_1*ER2;
    lam_OLS_1 = [lam_OLS_1 lam1(2)];
    
    X = [ones(N,1) cov_g cov_h2];
    M_1  =  eye(N) -  X*inv(X'*X)*X';
    P_1 = inv(X'*X)*X';
    lam1 = P_1*ER2;
    lam_OLS_2 = [lam_OLS_2 lam1(2)];



%% double lasso  and single lasso 
     
    model_ds = DS_TSCV(Ri, gt, ht,1,100,100,10);
    tstat_ds = model_ds.lambdag_ds/model_ds.se_ds;
    lambda_ds= model_ds.lambdag_ds;
    lam_DBL  = [lam_DBL lambda_ds];
    t_DBL    = [t_DBL ; tstat_ds];
    
    reccoff_DBL(i,model_ds.select) = 1;
    reccoff_DBL1(i,model_ds.sel1) = 1;
    lambda_ss = model_ds.lambdag_ss;
    tstat_ss  = model_ds.lambdag_ss/model_ds.se_ss;
    lam_SL    = [lam_SL  lambda_ss];
    t_SL      = [t_SL ;tstat_ss];


% HBIC
     opts1.tau=1;
     opts1.J = 100;
     opts1.tau0 = 1;
     opts1.isCV = 0;
     opts1.outputmode = 2;
     opts1.mu1 = size(ht,1)/N;
     opts1.delta = 1e-6;
     opts1.ifintercept = 1;
     opts1.initial=[];
     model_DL0_hbic  = DESDAR_CV(Ri,gt,ht,opts1,10,200);
     lam_Dl0_hbic=[lam_Dl0_hbic model_DL0_hbic.lambdag_ds];
     tstat_Dl0_hbic = model_DL0_hbic.lambdag_ds./model_DL0_hbic.se_ds;
     tstat_Sl0_hbic = model_DL0_hbic.lambdag_ss./model_DL0_hbic.se_ss;
     lam_SL0_hbic=[lam_SL0_hbic model_DL0_hbic.lambdag_ss];
     t_DL0_hbic  =[t_DL0_hbic tstat_Dl0_hbic];
     t_SL0_hbic  =[t_SL0_hbic tstat_Sl0_hbic];
     reccoff_Dl0_hbic(i,model_DL0_hbic.select) = 1;
     reccoff_Dl01_hbic(i,model_DL0_hbic.sel1) = 1;
     reccoff_Dl02_hbic(i,model_DL0_hbic.sel2) = 1;
     reccoff_Dl03_hbic(i,model_DL0_hbic.sel3) = 1;  
    
%% CV 
     opts2.tau=1;
     opts2.J = 400;
     opts2.tau0 = 1;
     opts2.isCV = 1;
     opts2.outputmode = 1;
     opts2.mu1 = size(ht,1)/N;
     opts2.delta = 1e-6;
     opts2.ifintercept =1;
     opts2.initial='default';
     model_DL0  = DESDAR_CV(Ri,gt,ht,opts2,10,200);
     lam_Dl0=[lam_Dl0 model_DL0.lambdag_ds];
     tstat_Dl0 = model_DL0.lambdag_ds./model_DL0.se_ds;
     tstat_Sl0 = model_DL0.lambdag_ss./model_DL0.se_ss;
     lam_SL0=[lam_SL0 model_DL0.lambdag_ss];

     t_DL0  =[t_DL0 tstat_Dl0];
     t_SL0  =[t_SL0 tstat_Sl0];

     reccoff_Dl0(i,model_DL0.select) = 1;
     reccoff_Dl01(i,model_DL0.sel1) = 1;
     reccoff_Dl02(i,model_DL0.sel2) = 1;
     reccoff_Dl03(i,model_DL0.sel3) = 1;  
%% 


end
