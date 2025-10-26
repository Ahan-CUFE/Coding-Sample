% This is the program to provide the simulation figure in the appendix and
% an exmaple for the asymptotic performance

% Part I: CLT plots
% Part II: exmaple for the asymptotic performance


clear; clc;
close all;

% choose p, n an T for the plots

T_values = [240, 360, 480];   % T的取值
N_values = [200, 300, 400];   % N的取值
P_values = [60, 80, 100];     % P的取值
Results_res = [];  
Results_RMSE = [];  
factor_type = 3;
for T = T_values
    for N = N_values
        row_bias = [];
        row_RMSE = [];
        for P = P_values
        simName = strcat('T_',num2str(T),'_n_',num2str(N),'_p_',num2str(P));
        load(strcat(simName))
        bias_DBL = mean(lam_DBL') - lambdag';
        RMSE_DBL= sqrt(mean((lam_DBL-lambdag)'.^2));
        
        bias_Dl0 = mean(lam_Dl0') - lambdag';
        RMSE_Dl0 = sqrt(mean((lam_Dl0-lambdag)'.^2));

        bias_Dl0_hbic = mean(lam_Dl0_hbic') - lambdag';
        RMSE_Dl0_hbic = sqrt(mean((lam_Dl0_hbic-lambdag)'.^2));

        %%

% 重新排序后的数据行
        new_row_res = [...
        P, T, N, ...
        bias_DBL(factor_type), bias_Dl0(factor_type),bias_Dl0_hbic(factor_type), ...
   ];

        new_row_RMSE = [...
        P, T, N, ...
         RMSE_DBL(factor_type), RMSE_Dl0(factor_type),RMSE_Dl0_hbic(factor_type), ...

        ];


        row_bias = [row_bias new_row_res];
        row_RMSE = [row_RMSE new_row_RMSE];
        %%
        % new_row = {P,T, N, bias_DBL(1), RMSE_DBL(1),bias_Dl0(1),RMSE_Dl0(1),bias_DBL(2), RMSE_DBL(2),bias_Dl0(2),RMSE_Dl0(2),bias_DBL(3), RMSE_DBL(3),bias_Dl0(3),RMSE_Dl0(3)};
      
        end
        Results_res = [Results_res; row_bias];  
        Results_RMSE= [Results_RMSE ; row_RMSE];  
    end
end
%%

Results_res = array2table(Results_res);
Results_RMSE= array2table(Results_RMSE);

Results_res.Properties.VariableNames = {...
    'P', 'T', 'N', ...
    'bias_DBL_usefull', 'bias_Dl0_usefull', 'bias_sl0_usefull', ...
    'P2', 'T2', 'N2', ...
    'bias_DBL_usefull2', 'bias_Dl0_usefull2', 'bias_sl0_usefull2', ...
    'P3', 'T3', 'N3', ...
    'bias_DBL_usefull3', 'bias_Dl0_usefull3', 'bias_sl0_usefull3', ...
};

Results_RMSE.Properties.VariableNames={...
    'P', 'T', 'N', ...   
    'RMSE_DBL_usefull', 'RMSE_Dl0_usefull','RMSE_sl0_usefull', ...
     'P2', 'T2', 'N2', ...   
    'RMSE_DBL_usefull2', 'RMSE_Dl0_usefull2','RMSE_sl0_usefull2', ...
    'P3', 'T3', 'N3', ...   
    'RMSE_DBL_usefull3', 'RMSE_Dl0_usefull3','RMSE_sl0_usefull3', ...,
    };


%%
% Results.Properties.VariableNames = {'P','T', 'N',  'bias_DBL1',' RMSE_DBL1','bias_Dl01','RMSE_Dl01','bias_DBL2',' RMSE_DBL2','bias_Dl02','RMSE_Dl02','bias_DBL3',' RMSE_DBL3','bias_Dl03','RMSE_Dl03'};


output_file = 'res_results_useless.csv';
writetable(Results_res, output_file);


output_file = 'rmse_results_useless.csv';
writetable(Results_RMSE, output_file);


% 
% 
% %% Part I: CLT plots
% 
% nbin = 50; lw=1;
% binCtrs = linspace(-6,6,nbin);
% binWidth=binCtrs(2)-binCtrs(1);
% 
% fig11 = figure;
% 
% % double selection
% 
% subplot(3,2,1)
% 
% counts=hist(lambdag_ds_std(:,1),binCtrs);
% prob = counts / (K * binWidth);
% h2=bar(binCtrs,prob,'hist');
% 
% set(h2,'FaceColor',[.6 .6 .6]);
% set(h2,'EdgeColor',[.6 .6 .6]);
% 
% xgrid = linspace(-4,4,81);
% pdfReal=pdf('Normal',-4:0.1:4,0,1);
% line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
% xlim([-5,5]);
% ylim([0,0.5]);
% title('DS: Uselful','FontSize',11)
% 
% subplot(3,2,3)
% counts=hist(lambdag_ds_std(:,2),binCtrs);
% prob = counts / (K * binWidth);
% h2=bar(binCtrs,prob,'hist');
% 
% set(h2,'FaceColor',[.6 .6 .6]);
% set(h2,'EdgeColor',[.6 .6 .6]);
% 
% xgrid = linspace(-4,4,81);
% pdfReal=pdf('Normal',-4:0.1:4,0,1);
% line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
% xlim([-5,5]);
% ylim([0,0.5]);
% title('DS: Redundant','FontSize',11)
% 
% subplot(3,2,5)
% counts=hist(lambdag_ds_std(:,3),binCtrs);
% prob = counts / (K * binWidth);
% h2=bar(binCtrs,prob,'hist');
% 
% set(h2,'FaceColor',[.6 .6 .6]);
% set(h2,'EdgeColor',[.6 .6 .6]);
% 
% xgrid = linspace(-4,4,81);
% pdfReal=pdf('Normal',-4:0.1:4,0,1);
% line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
% xlim([-5,5]);
% ylim([0,0.5]);
% title('DS: Useless','FontSize',11)
% 
% % single selection 
% 
% subplot(3,2,2)
% 
% counts=hist(lambdag_ss_std(:,1),binCtrs);
% prob = counts / (K * binWidth);
% h2=bar(binCtrs,prob,'hist');
% 
% set(h2,'FaceColor',[.6 .6 .6]);
% set(h2,'EdgeColor',[.6 .6 .6]);
% 
% xgrid = linspace(-4,4,81);
% pdfReal=pdf('Normal',-4:0.1:4,0,1);
% line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
% xlim([-5,5]);
% ylim([0,0.5]);
% title('SS: Uselful','FontSize',11)
% 
% subplot(3,2,4)
% counts=hist(lambdag_ss_std(:,2),binCtrs);
% prob = counts / (K * binWidth);
% h2=bar(binCtrs,prob,'hist');
% 
% set(h2,'FaceColor',[.6 .6 .6]);
% set(h2,'EdgeColor',[.6 .6 .6]);
% 
% xgrid = linspace(-4,4,81);
% pdfReal=pdf('Normal',-4:0.1:4,0,1);
% line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
% xlim([-5,5]);
% ylim([0,0.5]);
% title('SS: Redundant','FontSize',11)
% 
% subplot(3,2,6)
% counts=hist(lambdag_ss_std(:,3),binCtrs);
% prob = counts / (K * binWidth);
% h2=bar(binCtrs,prob,'hist');
% 
% set(h2,'FaceColor',[.6 .6 .6]);
% set(h2,'EdgeColor',[.6 .6 .6]);
% 
% xgrid = linspace(-4,4,81);
% pdfReal=pdf('Normal',-4:0.1:4,0,1);
% line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
% xlim([-5,5]);
% ylim([0,0.5]);
% title('SS: Useless','FontSize',11)
% 
% saveas(fig11,'Figure11','epsc')
% 
% %% Part II: exmaple for the asymptotic performance
% 
% % true values for lambdas
% lambdag = [lambdag(1);0;0]; 
% 
% % MC bias for symbol
% disp('MC bias')
% disp({'useful','redundant','useless'})
% disp(mean(lambdag_ds)-lambdag')
% disp(mean(lambdag_ss)-lambdag')
% disp(mean(lambdag_ns)-lambdag')
% 
% 
% % MC RMSE for symbol
% disp('MC RMSE')
% disp({'useful','redundant','useless'})
% disp(sqrt(mean((lambdag_ds-ones(K,1)*lambdag').^2)))
% disp(sqrt(mean((lambdag_ss-ones(K,1)*lambdag').^2)))
% disp(sqrt(mean((lambdag_ns-ones(K,1)*lambdag').^2)))
% 
% 
% 
