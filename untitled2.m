clear; clc; close all;

set(groot, 'defaultTextInterpreter', 'tex', ...
           'defaultAxesTickLabelInterpreter', 'tex', ...
           'defaultLegendInterpreter', 'tex');

%% ========== 1. 参数网格设置 ==========
% 参数范围（可根据需要调整）
K_min = 10;    % 开环增益最小值
K_max = 300;   % 开环增益最大值
a_min = 0.01;     % PD微分时间常数最小值（a=0为无PD的原系统）
a_max = 1;  % PD微分时间常数最大值
grid_num = 200; % 网格密度，数值越大越精细，计算越慢

% 生成二维参数网格
K_grid = linspace(K_min, K_max, grid_num);
a_grid = linspace(a_min, a_max, grid_num);
[K_mat, a_mat] = meshgrid(K_grid, a_grid);

% 预分配性能指标矩阵
sigma_map = zeros(size(K_mat));  % 超调量，单位%
ts_map    = zeros(size(K_mat));  % 调节时间（2%准则），单位s

%% ========== 2. 遍历参数网格，计算性能指标 ==========
for i = 1:numel(K_mat)
    K_val = K_mat(i);
    a_val = a_mat(i);
    
    % 构造开环传递函数：PD控制器串联原对象
    % 分子：K*(s+20)*(1+a*s) = K*(a*s² + (20a+1)s + 20)
    num = K_val * [a_val, 20*a_val + 1, 20];
    % 分母：s*(s²+24s+104) = s³ + 24s² + 104s
    den = [1, 24, 104, 0];
    
    G_open = tf(num, den);
    G_close = feedback(G_open, 1);  % 单位负反馈闭环
    
    % 计算阶跃响应性能指标（默认2%准则的调节时间）
    info = stepinfo(G_close);
    sigma_map(i) = info.Overshoot;
    ts_map(i)    = info.SettlingTime;
end

%% ========== 3. 绘制双热力图 ==========
figure('Color','w','Position',[100,100,1150,460]);

% --- 子图1：超调量热力图 ---
subplot(1,2,1);
imagesc(K_grid, a_grid, sigma_map);
colormap jet;
hold on;
% 标注超调=10%的临界线（题目要求≤10%）
contour(K_grid, a_grid, sigma_map, [10 10], ...
    'r--', 'LineWidth',1.6, 'DisplayName','\sigma% = 10%');

xlabel('开环增益 K');
ylabel('PD常数 Td');
title('闭环系统超调量分布');
axis xy;  % 纵轴从下到上递增
grid on;
legend('Location','northeast');
c = colorbar;
c.Label.String = '超调量 (%)';

% --- 子图2：调节时间热力图 ---
subplot(1,2,2);
imagesc(K_grid, a_grid, ts_map);
colormap jet;
hold on;
% 标注调节时间=1s的临界线（题目要求≤1s）
contour(K_grid, a_grid, ts_map, [1 1], ...
    'r--', 'LineWidth',1.6, 'DisplayName','t_s = 1s');

xlabel('开环增益 K');
ylabel('PD常数 Td');
title('闭环系统调节时间分布');
axis xy;
grid on;
legend('Location','northeast');
c = colorbar;
c.Label.String = '调节时间 (s)';

%% ========== 4. 查找满足性能要求的最优参数 ==========
% 要求：超调≤10% 且 调节时间≤1s
valid_mask = (sigma_map <= 10) & (ts_map <= 1);

if any(valid_mask(:))
    % 在可行域内找调节时间最短的参数点
    ts_valid = ts_map;
    ts_valid(~valid_mask) = inf;
    [min_ts, idx_opt] = min(ts_valid, [], 'all', 'linear');
    
    K_opt  = K_mat(idx_opt);
    a_opt  = a_mat(idx_opt);
    sig_opt = sigma_map(idx_opt);
    
    % 在两张图上标注最优参数点
    subplot(1,2,1);
    plot(K_opt, a_opt, 'r*', 'MarkerSize',14, 'LineWidth',1.8, 'DisplayName','最优参数点');
    subplot(1,2,2);
    plot(K_opt, a_opt, 'r*', 'MarkerSize',14, 'LineWidth',1.8, 'DisplayName','最优参数点');
    
    % 控制台输出结果
    fprintf('===== 同时满足 σ%%≤10%% 且 t_s≤1s 的最优参数（调节时间最短）=====\n');
    fprintf('开环增益 K  = %.2f\n', K_opt);
    fprintf('微分时间 a  = %.4f s\n', a_opt);
    fprintf('对应超调量  = %.2f %%\n', sig_opt);
    fprintf('对应调节时间 = %.4f s\n', min_ts);
else
    disp('当前参数范围内，不存在同时满足超调≤10%且调节时间≤1s的参数组合。');
end

%% ========== 5. 对比：无PD时原系统的性能 ==========
K_ref = 120;  % 可替换为第(4)问你选定的K值
num_ref = K_ref * [0, 1, 20];
den_ref = [1, 24, 104, 0];
G_ref = feedback(tf(num_ref, den_ref), 1);
info_ref = stepinfo(G_ref);

fprintf('\n===== 无PD控制器（a=0）时，K=%.2f 的原系统性能 =====\n', K_ref);
fprintf('超调量   = %.2f %%\n', info_ref.Overshoot);
fprintf('调节时间 = %.4f s\n', info_ref.SettlingTime);
