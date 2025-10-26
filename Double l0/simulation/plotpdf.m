% data = (lam_Dl0(1,:)-lambdag(1))./std(lam_Dl0(1,:)-lambdag(1));
% data = (lam_DBL(1,:)-lambdag(1))./std(lam_DBL(1,:)-lambdag(1));
% data = (lam_SL0(1,:)-lambdag(1))./std(lam_SL0(1,:)-lambdag(1));
% histogram(data, 'Normalization','pdf',...
% 'FaceColor',[0.7 0.7 0.7],'EdgeColor','none'); % 灰色填充、无边框
% hold on
% x = -5:0.01:5;
% plot(x, normpdf(x,0,1), 'k--','LineWidth',1.5) % 标准正态曲线
% title('DS: Useful')


%% Part I: CLT plots
lw = 1 ;
fig11 = figure;
% double selection
subplot(3,2,1)
% load('t_bias_rec.mat');
std_Dl0_usefull = t_bias_Dl0(1,:);


% std_Dl0_usefull = (lam_Dl0(1,:)-lambdag(1))./std(lam_Dl0(1,:)-lambdag(1));
histogram(std_Dl0_usefull, 'Normalization','pdf',...
'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','BinWidth',0.2); % 灰色填充、无边框


xgrid = linspace(-4,4,81);
pdfReal=pdf('Normal',-4:0.1:4,0,1);
line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
xlim([-5,5]);
ylim([0,0.5]);
title('Dl0: Uselful','FontSize',11)




subplot(3,2,3)

% std_Dl0_red = (lam_Dl0(2,:)-lambdag(2))./std(lam_Dl0(2,:)-lambdag(2));
std_Dl0_red = t_bias_Dl0(2,:);
histogram(std_Dl0_red, 'Normalization','pdf',...
'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','BinWidth',0.2); % 灰色填充、无边框
xgrid = linspace(-4,4,81);
pdfReal=pdf('Normal',-4:0.1:4,0,1);
line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
xlim([-5,5]);
ylim([0,0.5]);
title('Dl0: Redundant','FontSize',11)




subplot(3,2,5)

% std_Dl0_useless = (lam_Dl0(3,:)-lambdag(3))./std(lam_Dl0(3,:)-lambdag(3));
std_Dl0_useless = t_bias_Dl0(3,:);
histogram(std_Dl0_useless, 'Normalization','pdf',...
'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','BinWidth',0.5); % 灰色填充、无边框
xgrid = linspace(-4,4,81);
pdfReal=pdf('Normal',-4:0.1:4,0,1);
line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
xlim([-5,5]);
ylim([0,0.5]);
title('Dl0: Useless','FontSize',11)

%%
% single selection 

subplot(3,2,2)

% std_sl0_usefull = (lam_SL0(1,:)-lambdag(1))./std(lam_SL0(1,:)-lambdag(1));
std_sl0_usefull = t_bias_Dl0_hbic(1,:);
histogram(std_sl0_usefull, 'Normalization','pdf',...
'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','BinWidth',0.2); % 灰色填充、无边框

xgrid = linspace(-4,4,81);
pdfReal=pdf('Normal',-4:0.1:4,0,1);
line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
xlim([-5,5]);
ylim([0,0.5]);
title('Dl0_{HBIC}: Uselful','FontSize',11)




subplot(3,2,4)

% std_sl0_Redundant = (lam_SL0(2,:)-lambdag(2))./std(lam_SL0(2,:)-lambdag(2));
std_sl0_Redundant = t_bias_Dl0_hbic(2,:);
histogram(std_sl0_Redundant, 'Normalization','pdf',...
'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','BinWidth',0.2); % 灰色填充、无边框

xgrid = linspace(-4,4,81);
pdfReal=pdf('Normal',-4:0.1:4,0,1);
line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
xlim([-5,5]);
ylim([0,0.5]);

title('Dl0_{HBIC}: Redundant','FontSize',11)






subplot(3,2,6)
% std_sl0_Useless = (lam_SL0(3,:)-lambdag(3))./std(lam_SL0(3,:)-lambdag(3));
std_sl0_Useless = t_bias_Dl0_hbic(3,:);
histogram(std_sl0_Useless, 'Normalization','pdf',...
'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','BinWidth',0.2); % 灰色填充、无边框

xgrid = linspace(-4,4,81);
pdfReal=pdf('Normal',-4:0.1:4,0,1);
line(xgrid,pdfReal,'color','k','linestyle','--','linewidth',lw);
xlim([-5,5]);
ylim([0,0.5]);

title('Dl0_{HBIC}: Useless','FontSize',11)


%%



saveas(fig11,'Figure11','png')




