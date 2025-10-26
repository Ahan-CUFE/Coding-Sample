% % 假设 reccoff_DBL、reccoff_Dl0、reccoff_Dl0_hbic 都是 [1000 × p]
% [numSim, p] = size(reccoff_DBL);
% 
% % 统计每个因子被选中的次数
% count_DBL      = sum(reccoff_DBL     ~= 0, 1);   % DBL selection
% count_Dl0      = sum(reccoff_Dl0     ~= 0, 1);   % Dl0 selection
% count_DPL      = sum(reccoff_DBL_hbic~= 0, 1);
% count_Dl0_hbic = sum(reccoff_Dl0_hbic~= 0, 1);   % Dl0 (HBIC)
% [numSim, p] = deal(1000, length(count_DBL));  
% [numSim, p] = deal(1000, length(count_DBL));  
% 
% x = 1:p;
% grayColor = [0.7 0.7 0.7];
% barLW     = 0.5;    % 柱子边框宽度
% lineLW    = 1.5;    % 区分线宽度
% sepPos    = 4.5;    % 在第 4 和第 5 因子之间
% 
% figure;
% 
% %% 第一张：DBL
% subplot(3,1,1);
% bar(x, count_DBL, 'FaceColor', grayColor, 'EdgeColor','none', 'LineWidth',barLW);
% hold on;
% % 划一条红色虚线在 x = 4.5 处
% xline(sepPos, '--r', 'LineWidth', lineLW);
% hold off;
% xlim([0 p+1]);   ylim([0 numSim]);
% ylabel('Selection Count');
% title('1^{st} Selection: DBL','FontSize',11);
% set(gca,'FontSize',9);
% 
% %% 第二张：Dl0
% subplot(3,1,2);
% bar(x, count_Dl0, 'FaceColor', grayColor, 'EdgeColor','none', 'LineWidth',barLW);
% hold on;
% xline(sepPos, '--r', 'LineWidth', lineLW);
% hold off;
% xlim([0 p+1]);   ylim([0 numSim]);
% ylabel('Selection Count');
% title('2^{nd} Selection: Dl0','FontSize',11);
% set(gca,'FontSize',9);
% 
% %% 第三张：Dl0_{HBIC}
% subplot(3,1,3);
% bar(x, count_Dl0_hbic, 'FaceColor', grayColor, 'EdgeColor','none', 'LineWidth',barLW);
% hold on;
% xline(sepPos, '--r', 'LineWidth', lineLW);
% hold off;
% xlim([0 p+1]);   ylim([0 numSim]);
% ylabel('Selection Count');
% xlabel('Factor ID (simulated)');
% title('Total Selection: Dl0_{HBIC}','FontSize',11,'Interpreter','tex');
% set(gca,'FontSize',9);
% 
% % （可选）微调子图间距
% set(gcf,'Position',[100 100 600 700]);


% 假设 reccoff_DBL、reccoff_Dl0、reccoff_DBL_hbic、reccoff_Dl0_hbic 都是 [1000 × p]
[numSim, p] = size(reccoff_DBL);

% 统计每个因子被选中的次数
count_DBL      = sum(reccoff_DBL      ~= 0, 1);   % DBL selection
count_Dl0      = sum(reccoff_Dl0      ~= 0, 1);   % Dl0 selection
count_DPL      = sum(reccoff_DBL_hbic ~= 0, 1);   % Double Post-Lasso (DBL_HBIC)
count_Dl0_hbic = sum(reccoff_Dl0_hbic ~= 0, 1);   % Dl0 (HBIC)

x       = 1:p;
grayColor = [0.7 0.7 0.7];
barLW   = 0.5;    % 柱子边框宽度
lineLW  = 1.5;    % 区分线宽度
sepPos  = 4.5;    % 在第 4 和第 5 因子之间

figure('Position',[100 100 600 900]);

%% 第一张：DBL
subplot(4,1,1);
bar(x, count_DBL, 'FaceColor', grayColor, 'EdgeColor','none','LineWidth',barLW);
hold on;
xline(sepPos, '--r', 'LineWidth', lineLW);
hold off;
xlim([0 p+1]);   ylim([0 numSim]);
ylabel('Selection Count');
title('1^{st} Selection: Double Lasso','FontSize',11);
set(gca,'FontSize',9);

%% 第二张：Dl0
subplot(4,1,2);
bar(x, count_Dl0, 'FaceColor', grayColor, 'EdgeColor','none','LineWidth',barLW);
hold on;
xline(sepPos, '--r', 'LineWidth', lineLW);
hold off;
xlim([0 p+1]);   ylim([0 numSim]);
ylabel('Selection Count');
title('2^{nd} Selection: Double SDAR','FontSize',11);
set(gca,'FontSize',9);

%% 第三张：Double Post-Lasso
subplot(4,1,3);
bar(x, count_DPL, 'FaceColor', grayColor, 'EdgeColor','none','LineWidth',barLW);
hold on;
xline(sepPos, '--r', 'LineWidth', lineLW);
hold off;
xlim([0 p+1]);   ylim([0 numSim]);
ylabel('Selection Count');
title('3^{rd} Selection: Double Post‐Lasso with HBIC','FontSize',11,'Interpreter','tex');
set(gca,'FontSize',9);

%% 第四张：Dl0_{HBIC}
subplot(4,1,4);

bar(x, count_Dl0_hbic, 'FaceColor', grayColor, 'EdgeColor','none','LineWidth',barLW);
hold on;
xline(sepPos, '--r', 'LineWidth', lineLW);
hold off;
xlim([0 p+1]);   ylim([0 numSim]);
ylabel('Selection Count');
xlabel('Factor ID (simulated)');
title('4^{th} Selection: Double SDAR with HBIC','FontSize',11,'Interpreter','tex');
set(gca,'FontSize',9);


% （可选）微调子图间距
set(gcf,'Position',[100 100 600 900]);