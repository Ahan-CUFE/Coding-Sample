
% % factor
% allfactors = csvread('factors.csv',1,0);
% date = allfactors(:,1);
% rf = allfactors(:,2);
% factors = allfactors(:,3:end);
% 
% L = length(date);
% P = size(factors,2);
% 
% % test portfolios
% port_3x2 = csvread('port_3x2.csv',0,1);
% port_3x2 = port_3x2 - rf*ones(1,size(port_3x2,2)); % excess return
% 
% % other information
% summary = readtable('summary.csv');
% factorname = summary.Row;
% factorname_full = summary.Descpription;
% year_pub = summary.Year;
% year_end = summary.Year_end;
% 
% port_3x2_id = readtable('port_3x2_id.csv');
% 
% mkt_ind = find(ismember(factorname,'MktRf'));
% smb_ind = find(ismember(factorname,'SMB'));
% hml_ind = find(ismember(factorname,'HML'));
% rmw_ind = find(ismember(factorname,'RMW'));
% cma_ind = find(ismember(factorname,'CMA'));
% kk = 10; % minimun number of stocks in a portfolio
% 
% include_3x2 = find(port_3x2_id.min_stk6>=kk)';
% port_3x2b = [];
% 
% for i = 1:P
%     if ismember(i,include_3x2)
%         port_3x2b = [port_3x2b port_3x2(:,(i*6-5):(i*6))];
%     end
% end
% 
% Ri = port_3x2b; % test asset
% % choose control factors before 2012
% ContrlList = find(year_pub < 2012);
% ControlFactor = factors(:,ContrlList);
% FF3 = factors(:,[mkt_ind,smb_ind,hml_ind])';
% FF5 = factors(:,[mkt_ind,smb_ind,hml_ind,rmw_ind,cma_ind])';
% % test factors since 2012
% TestList = find(year_pub >= 2012);
% TestFactor = factors(:,TestList);
% 
% 
% Useless_index = setdiff([1:150], [mkt_ind,smb_ind,hml_ind,rmw_ind,cma_ind])
% Uselessfactor = factors(:,Useless_index)'
% 
% save('/Users/gaoyihan/Desktop/factor_zoo/code/doublel0/para_spefic.mat','Uselessfactor','FF5','Ri')
load('para_spefic.mat')
%%
% ht1 = FF5([1,3,4,5],:);
% gt = FF5(2,:);  % donate SMB as usefull one .

ht1 = FF5([1,2,3,4],:);
gt = FF5(5,:);

[m, n] = size(ht1');
X = ht1';
for i = 1:n
    for j = 1:i-1
        % 将第i列变量与前面所有列正交
        X(:,i) = X(:,i) - (X(:,i)'*X(:,j)) * X(:,j) / (X(:,j)'*X(:,j));
    end
end
ht1 = X';

%%
% gt = gt - mean(gt);
% ht1 = ht1 - mean(ht1')';
% ht=FF5;
ht2=Uselessfactor(1:132,:);
%%
%
sigma2_target = 0.001;

% 创建一个相同尺寸的矩阵用于存储处理后的数据
X = ht2;
X_new = zeros(size(X));

% nomissing = (sum(isnan(X_new),1)==0)';
% 对每个特征（每一行）进行处理
for i = 1:size(X,1)
    % 提取第 i 个特征的数据行
    Xi = X(i,:); 
    
    % 计算该特征的均值和方差
    mu = nanmean(Xi);
    sigma2 = nanvar(Xi);
    
    % 若原方差为零(或极小)，可能需要特殊处理，以避免除0错误
    if sigma2 == 0
        % 若方差为0，数据全相等，将此行数据略微扰动或直接赋值为常量
        % 这里假设原数据为常量，全为 mu
        % 可直接赋值为相同的mu（方差依然为0）或添加微小扰动
        X_new(i,:) = mu; % 无法创建0.001方差，除非添加噪声
    else
        % 根据目标方差计算缩放因子
        a = sqrt(sigma2_target / sigma2);

        % 对此特征进行缩放变换
        % 确保变换后均值不变，方差为目标值
        X_new(i,:) = mu + a*(Xi - mu);
    end
end

% 验证处理结果：检查每一行的新方差
for i = 1:size(X_new,1)
    new_var = var(X_new(i,:));
    fprintf('Feature %d: new variance = %.5f\n', i, new_var);
end

%%


ht2 = X_new;



ht=[ht1;ht2];
% gt =Uselessfactor(133:145,:);
Ri = Ri';
% data information
n = size(Ri,1);
p1 = size(ht1,1);
p2 = size(ht2,1);
dg = size(gt,1);


tmp1  = nancov([ht2',Ri']);
Ch2 = tmp1((p2+1):end,1:p2);%useless

tmp2  = nancov([ht1',Ri']);
Ch1 = tmp2((p1+1):end,1:p1);%usefull

tmp3  = nancov([gt',Ri']);
Cg = tmp3((dg+1):end,1:dg);

ER    = nanmean(Ri,2);


%% normliza the level 
% p=p1;
% beta = NaN(n,p);
% for i = 1:p
%   beta(:,i) = Ch1(:,i)/nanvar(ht1(i,:));
% end
% penalty = nanmean(beta.^2,1);
% penalty = penalty./nanmean(penalty); % normalize the level
% 
% Ch1 = Ch1*(diag(penalty));
% 
% p=p2;
% beta = NaN(n,p);
% for i = 1:p
%   beta(:,i) = Ch2(:,i)/nanvar(ht2(i,:));
% end
% penalty = nanmean(beta.^2,1);
% penalty = penalty./nanmean(penalty); % normalize the level
% 
% Ch2 = Ch2*(diag(penalty));
% %如此矫正可能会因为换手率与beta的相关性过高导致数值计算溢出
%%



%ch1

mean_ch1 = nanmean(Ch1);
sigma_ch1 = cov(Ch1);

cov_h10 =[ones(n,1) Ch1];
M_Ch0 = eye(n) - cov_h10*inv(cov_h10'*cov_h10)*cov_h10';
P_Ch0 = inv(cov_h10'*cov_h10)*cov_h10';
%C_eps
C_eps = M_Ch0*Ch2;
sigma_C_eps = cov(C_eps);
mean_eps = nanmean(C_eps);
%theta
lambda_ch0 = P_Ch0*Ch2;
theta0 = lambda_ch0(1,:);
theta1 = lambda_ch0(2:end,:);

%theta0:
mean_theta0=nanmean(theta0);
sigma_theta0 = cov(theta0);

mean_theta1 = nanmean(theta1');
sigma_theta1= cov(theta1');
%%
%only simulate the cov of Ce and the level of chi,then we will construt chi
%randomly with different information strength.


Ch = [Ch1 Ch2];
Ch0 = [ones(n,1) Ch];

Ch10=[ones(n,1) Ch1];

M_ch01 = eye(n) - Ch10*inv(Ch10'*Ch10)*Ch10';
P_Ch0 = inv(Ch10'*Ch10)*Ch10';

Ce1 = M_ch01*Cg;                   %求出Ce1，Ce2自由生成
chi0 = P_Ch0*Cg;
xi = chi0(1,:);

mean_Ce =nanmean(Ce1);
sigma_ce = cov(Ce1);


mean_xi = nanmean(xi);%xi 代表两个真实值均值的差异，不可以改变。故xi不用生成
% sigma_xi = cov(xi);

chi = chi0(2:end,:);
chi1 = chi(1:end,:);%to give the usefull factor' loading to cov_g
% chi2 = chi(6:end,:);%to give the useless factor' loading to cov_g ,chi2
% shou given by persion.

mean_chi1 = nanmean(nanmean(chi1,2));%chi1 and chi2 实际上应是固定的，不能有随机性
sigma_chi1 = cov(nanmean(chi1,2));
% nanmean_chi2 = nanmean(nanmean(chi2,2));
% sigma_chi2 = cov(nanmean(chi2,2));


%Cz sigma_Cz eta 
P_h=ht1'*inv(ht1*ht1');
% M_ch=eye(n)-Ch1*inv(Ch1'*Ch1)*Ch1';
% M_h=eye(n)-ht1*inv(ht1'*ht1)*ht1';

eta = gt*P_h; %to give eta
% eta1 = eta(1:4,:);
eta1 = zeros(p1+p2,1);
eta1(1:p1) = eta;

% eta2 = eta(6:end,:);
% nanmean_eta1 = nanmean(nanmean(eta1));  %eta here is etaT in support
% nanmean_eta2 = nanmean(nanmean(eta2));
% sigma_eta1 = cov(nanmean(eta1,2));
% sigma_eta2 = cov(nanmean(eta2,2));

Cz = Cg - Ch1*eta';
zt = gt - eta*ht1;
sigma_z = cov(zt);%%%%
mean_z = nanmean(zt);
%gamma0 lambdag lambdah 

CgCh0 = [ones(n,1) Cg Ch1];


P_CgCh0 = inv(CgCh0'*CgCh0)*CgCh0';


% to generate non-zero lambda_g and lambda_h
lambda_all = P_CgCh0*ER;
lambda0 = lambda_all(1);
lambda_g = nanmean(lambda_all(2));
lambda_h1 = lambda_all(3:6);
% lambda_h = nanmean(lambda_all(3:end));
sigma_lambda_h1 = cov(lambda_all(3:6));%lambdah1不可随机！


%calculate ut
mean_h = nanmean([ht1 ;ht2]');
sigma_h = nancov([ht1 ;ht2]');
% diag1 = 1.5e-5*eye(size(sigma_h));
% sigma_h = sigma_h + diag1;
sigma_h1 = cov(ht1');




ratio_sigma_h2 = 0.001;
sigma_h2 = ratio_sigma_h2 * diag(ones(p2,1));

betag = Cz*inv(sigma_z);

betah = Ch1*inv(sigma_h1)-betag*eta;
ut = Ri - ER -betag*gt-betah*ht1;
sigma_ut = nancov(mean(ut));

% 
save('/Users/gaoyihan/Desktop/factor_zoo/code/DBL0/doublel0/simulation2/simudata/hyparamter2.mat','mean_eps','sigma_C_eps','mean_theta0','sigma_theta0', ...
 'mean_theta1','sigma_theta1' ,'sigma_C_eps' ,'mean_eps','mean_Ce','sigma_ce','xi','mean_xi', 'theta1','theta0'...
,'mean_chi1','sigma_chi1','mean_z', ...
'sigma_z','lambda0','lambda_g','lambda_h1','sigma_lambda_h1','sigma_ut','mean_ch1','sigma_ch1','sigma_h','sigma_h1','sigma_h2','mean_h','chi1','eta1','eta')

% 
% clc;
% clear;


