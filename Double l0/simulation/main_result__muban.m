% 参数初始化
T_values = [240, 360, 480];   % T的取值
N_values = [200, 300, 400];   % N的取值
P_values = [60, 80, 100];     % P的取值

% 初始化结果存储表
Results = table();  

% 文件名模板
file_template = 'new_T_%d_n_%d_p_%d.mat';  % 文件名格式

% 循环读取文件并计算误差
for T = T_values
    for N = N_values
        for P = P_values
            % 构建文件名
            file_name = sprintf(file_template, T, N, P);
            if exist(file_name, 'file')  % 检查文件是否存在
                % 加载文件
                data = load(file_name);
                
                % 检查是否包含必要变量
                if isfield(data, 'lam_DBL') && isfield(data, 'lam_Dl0') && isfield(data, 'lambdag')
                    % 计算误差
                    DBL_error = mean(data.lam_DBL') - data.lambdag';
                    Dl0_error = mean(data.lam_Dl0') - data.lambdag';
                    
                    % 将结果追加到表中
                    new_row = {T, N, P, DBL_error, Dl0_error};
                    Results = [Results; new_row];  
                else
                    fprintf('警告: 文件 %s 中缺少必要变量.\n', file_name);
                end
            else
                fprintf('警告: 文件 %s 未找到.\n', file_name);
            end
        end
    end
end

% 设置表头
Results.Properties.VariableNames = {'T', 'N', 'P', 'DBL_Error', 'Dl0_Error'};

% 导出表格到CSV文件
output_file = 'error_results.csv';
writetable(Results, output_file);

fprintf('所有误差计算完成，结果已保存到文件: %s\n', output_file);