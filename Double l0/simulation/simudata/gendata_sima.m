function [covsigmaz,sigmah,eta,sigma_ut,betag,betah,ER,lambdag,chi1,betag1,betah1 ]=gendata_sima(N,T,d,p1,p2)
%========================================================================
% Input:
% N : num of asset
% T : len of time
% d : num of factor of  interest
% p : num of control factor
% Output:
% Ri: asset price
% gt: factor interest
% ht: control factor

%=======================================================================
rng(100);
load('hyparamter2.mat');
p=p1+p2;
n=N;
Ri = nan(N,T);
gt = nan(d,T);
ht = nan(p,T);

%because we can test factor indivadully, so it is nature to just consider
%one usefull factor here with Ce calculated form FF5 , the Ce of other d-1 interesting
%factor we can simulation form mean 0 norm distribution with noise
%level like ce

Ce_1 = mvnrnd(mean_Ce,sigma_ce,n);

Ce_2 = mvnrnd(zeros(1,d-1),diag(sigma_ce*ones(d-1,1)),n);
Ce = [Ce_1 Ce_2];

%here we boostrap uniformly




% if p1 <= size(sigma_ch1,1) % 有效控制元素小于Ch的维度时 ， 进行向下取样
%     index  = sort(randperm(size(sigma_ch1,1), p1));
%     Ch1 = mvnrnd(mean_ch1(index),sigma_ch1(index,index),n);
% end
% 
% %gen C_eps n*p2 dim
% if p2 <= size(sigma_C_eps,1) % 有效控制元素小于Ch的维度时 ， 进行向下取样
%     index  = sort(randperm(size(sigma_C_eps,1), p2));
%     C_eps = mvnrnd(mean_eps(index),sigma_C_eps(index,index),n);
% end

%%
if p1 <= size(sigma_ch1,1) % 有效控制元素小于Ch的维度时 ， 进行向下取样
    indexh1  = sort(randperm(size(sigma_ch1,1), p1));
    Ch1 = mvnrnd(mean_ch1(indexh1),sigma_ch1(indexh1,indexh1),n);
end

if p2 <= size(sigma_C_eps,1) % 有效控制元素小于Ch的维度时 ， 进行向下取样
    indexh2  = sort(randperm(size(sigma_C_eps,1), p2));
    C_eps = mvnrnd(mean_eps(indexh2),sigma_C_eps(indexh2,indexh2),n);
end



%%
%gen theta0 and theta1
theta0 = theta0([indexh2])';
theta1 = theta1(indexh1,indexh2)';
% theta0 = mvnrnd(mean_theta0,sigma_theta0,p2);
% if p1 <= size(sigma_theta1,1) % 有效控制元素小于Ch的维度时 ， 进行向下取样
%     % index  = sort(randperm(size(sigma_theta1,1), p1));
%     index = indexh1;
%     theta1 = mvnrnd(mean_theta1(index),sigma_theta1(index,index),p2);
% end
%calculate Ch2
Ch2 = ones(n,1)*theta0'+Ch1*theta1'+C_eps;


%%






%%







%calculate Cg
Ch = [Ch1 Ch2];
xi=[mean_xi,zeros(1,d-1)];
%% chi and eta
chi1 = chi1(indexh1);
% chi1 = mvnrnd(mean_chi1,sigma_chi1,p1);  % ht1 对 gt1 的解释  

chi2 = nan(p1,d-1);
% 
% for i =1 : d-1
%     chi2(:,i) = mvnrnd(mean_chi2,sigma_chi2 , p1);% ht1 对gt2的解释 4*（d-1） 之后还要转秩
% end
chi2 = [0.03,0,0.06,0;0,0.0,0,0.0]';


chi = zeros(d,p1+p2);
chi(1,1:p1) = chi1;
chi(2:d,1:p1) = chi2';
Cg = ones(n,1)*xi + Ch*chi'  +Ce;
%%

%% eta

% eta1 = mvnrnd(mean_eta1,sigma_eta1,p1);
eta1 = eta(indexh1);
eta2 = nan(p1,d-1);
% for i = 1:d-1
%     eta2(:,i) = mvnrnd(mean_eta2,sigma_eta2,p1);% ht1 对gt2的解释 4*（d-1） 之后还要转秩
% end

eta2 = chi2;

eta = zeros(d,p1+p2);
eta(1,1:p1) =eta1;
eta(2:d,1:p1) = eta2';

Cz = Cg - Ch*eta';






lambdag = zeros(d,1);
%generate information strength

lambdag(1) = lambda_g;
% lambdah1 = mvnrnd(lambda_h,sigma_lambda_h,p1);
lambdah1 = lambda_h1(indexh1);
lambdah = zeros(p1+p2,1);
lambdah(1:p1) = lambdah1; 
ER = ones(n,1)*lambda0 + Cg*lambdag + Ch*lambdah;


covsigmaz = diag(sigma_z*ones(d,1));
betag = Cz*inv(covsigmaz);

if (p1+p2)<=size(sigma_h,1) && p1<=size(sigma_h1,1)&&p2<=size(sigma_h2,1)
    % index_1 = sort(randperm(size(sigma_h1,1), p1));
    % index_2 = size(sigma_h1,1)+sort(randperm(size(sigma_h2,1), p2));
    index_1 = indexh1;
    index_2 = 4+indexh2;
    index = [index_1 index_2];
    sigmah = sigma_h(index,index);
end

betah = Ch*inv(sigmah)-betag*eta;


%%

%betag1 and betah1

Cz1 = Cg(:,1) - Ch1*eta1';

betag1 = Cz1*inv(sigma_z);
betah1 = Ch1*inv(sigmah(index_1,index_1)) - betag1*eta1;




%gen ht,zt and calculate rt 
% zt = nan(d,T);
% ut = nan(n,T);
% meanz = zeros(d,1);
% % meanz(1) = mean_z;
% for t = 1:T
%     ztht = mvnrnd(zeros(p1+p2+d,1),blkdiag(covsigmaz, sigmah),1);
%     % ht(:,t) = mvnrnd(zeros(p1+p2,1),sigmah,1);
%     % ht(:,t) = mvnrnd(mean_h([indexh1 4+indexh2]),sigmah,1);
%     % zt(:,t) = mvnrnd(zeros(d,1),covsigmaz,1);
%     zt(:,t) = ztht(1:d);
%     ht(:,t) = ztht(d+1:end);
%     % zt(:,t) = mvnrnd(meanz,covsigmaz,1);
%     gt(:,t) = eta*ht(:,t)+zt(:,t);
%     % gt(:,t) = eta*ht(:,t)+;
%     if n <=size(sigma_ut,1)
%         indexu  = sort(randperm(size(sigma_ut,1), n));
%         % ut(:,t) = mvtrnd(sigma_ut(indexu,indexu),4,1);
% 
%     end
% %     ut(:,t) = mvtrnd(sigma_ut*diag(ones(n,1)),5,1);
%     Ri(:,t) = ER +betag*gt(:,t) + betah*ht(:,t);%+ut(:,t);
% end
% 
% 
% save('/Users/gaoyihan/Desktop/factor_zoo/code/new_code/hyparamter_genedata.mat','covsigmaz','sigmah','eta','sigma_ut', ...
%    'betag','betah' )
% 
% 
% 
% ztht = mvnrnd(zeros(p1+p2+d,1),blkdiag(covsigmaz, sigmah),T);
% zt = ztht(:,1:d)';
% ht= ztht(:,d+1:end)';
% gt = eta*ht+zt;
% ut = mvtrnd(sigma_ut*diag(ones(n,1)),5,T);
% 
% Ri= ER +betag*gt + betah*ht+ut';
% 
% 
% 
% ER2 = mean(Ri,2);
% gt1 = gt(1,:);
% ht1 = ht(1:4,:);
% 
% tmp1 = nancov([gt1;Ri]');
% cov_g = tmp1(2:end,1);
% 
% tmp2 = nancov([ht1;Ri]');
% cov_h = tmp2(5:end,1:4);
% 
% X = [ones(N,1) cov_g cov_h];
% M_1  =  eye(N) -  X*inv(X'*X)*X';
% P_1 = inv(X'*X)*X';
% lam1 = P_1*ER2;
% lam11 = lam1(2);






