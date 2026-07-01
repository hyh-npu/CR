%% 绘制原开环传递函数根轨迹图
clear; close; clc;
s= tf('s');
gs = (s+20)/(s*(s*s + 24*s +104));
%rlocus(gs);
step(feedback(20.3*gs,1));

%% 根轨迹分析
function rlocus_annotated(sys)

% RLOCUS_ANNOTATED 绘制根轨迹并自动标注渐近线、分离点、虚轴交点、出射/入射角
% 输入参数：
%   sys - 开环传递函数（tf 对象，单位负反馈结构）
% 标注内容：
%   1. 根轨迹渐近线（虚线）+ 重心
%   2. 实轴分离/会合点
%   3. 根轨迹与虚轴的交点及对应增益
%   4. 复极点的出射角
%   5. 复零点的入射角
% 依赖：Symbolic Math Toolbox（用于解析求解关键点）

% 1. 基础参数提取
[num, den] = tfdata(sys, 'v');  
p = pole(sys);   % 全部开环极点
z = zero(sys);   % 全部开环零点
n = length(p);
m = length(z);

% 2. 绘制根轨迹基底
figure('Color','w');
rlocus(sys);
hold on;
axis equal;  
grid on;
title('根轨迹图');
xlabel('实轴 Re(s)');
ylabel('虚轴 Im(s)');

% ====================== 3. 绘制渐近线 ======================
if n > m
    n_m = n - m;
    % 渐近线重心
    sigma_a = (sum(p) - sum(z)) / n_m;
    % 渐近线夹角（角度制）
    phi_deg = ((0:n_m-1)*2 + 1) * 180 / n_m;
    phi_rad = deg2rad(phi_deg);
    
    % 自适应射线长度
    x_range = xlim;
    t_max = max(abs(x_range)) * 2.5;
    t = linspace(0, t_max, 300);
    
    % 绘制每条渐近线
    for k = 1:n_m
        s_asym = sigma_a + t * exp(1j * phi_rad(k));
        plot(real(s_asym), imag(s_asym), 'k--', ...
            'LineWidth', 1.2, 'HandleVisibility','off');
    end
    
    % 标注重心
    plot(real(sigma_a), imag(sigma_a), 'rs', ...
        'MarkerFaceColor','r', 'DisplayName','渐近线重心');
    text(real(sigma_a)+0.3, imag(sigma_a)+0.3, ...
        sprintf('σ_a=%.2f', real(sigma_a)), 'Color','r');
end

% ====================== 4. 计算并标注分离/会合点 ======================
syms s_sym
% 分离点方程
lhs = sum(1 ./ (s_sym - p));
rhs = sum(1 ./ (s_sym - z));
eq = lhs == rhs;
s_sep_all = double(solve(eq, s_sym));

% 筛选有效点：K>0 且为实根
sep_points = [];
sep_K = [];
for k = 1:length(s_sep_all)
    s_val = s_sep_all(k);
    if ~isfinite(s_val) || isnan(s_val)
        continue;
    end
    % 仅保留近似实数解
    if abs(imag(s_val)) < 1e-6
        s_real = real(s_val);
        % 计算对应的增益 K
        L_at_s = evalfr(sys, s_real);
        % 避免除以零或极小值
        if abs(L_at_s) < 1e-12
            continue;
        end
        K_val = -1 / L_at_s;
        K_val = real(K_val); %去掉数值误差的虚部
        if K_val > 1e-6
            sep_points(end+1) = s_real;
            sep_K(end+1) = K_val; 
        end
    end
end

if ~isempty(sep_points)
    plot(sep_points, zeros(size(sep_points)), 'mo', ...
        'MarkerFaceColor','m', 'DisplayName','分离/会合点');
    for i = 1:length(sep_points)
        text(sep_points(i), 0.5, ...
            sprintf('s=%.2f\nK=%.2f', sep_points(i), sep_K(i)), ...
            'Color','m', 'HorizontalAlignment','center');
    end
end
%% ====================== 5. 计算并标注虚轴交点 ======================
syms w_sym K_sym s_sym
s_jw = 1j * w_sym;
% 将系数向量转换为符号多项式，再在 s = jω 处代入
num_sym = poly2sym(num, s_sym);
den_sym = poly2sym(den, s_sym);
char_eq = subs(den_sym, s_sym, s_jw) + K_sym * subs(num_sym, s_sym, s_jw);
real_eq = real(char_eq) == 0;
imag_eq = imag(char_eq) == 0;


% 求解方程组
sol_jw = solve([real_eq, imag_eq], [w_sym, K_sym]);
w_vals = double(sol_jw.w_sym);
K_vals = double(sol_jw.K_sym);

% 筛选有效解：ω>0、K>0
valid_idx = (w_vals > 1e-6) & (K_vals > 1e-6);
w_cross = w_vals(valid_idx);
K_cross = K_vals(valid_idx);

if ~isempty(w_cross)
    plot(zeros(size(w_cross)),  w_cross, 'bo', ...
        'MarkerFaceColor','b', 'DisplayName','虚轴交点');
    plot(zeros(size(w_cross)), -w_cross, 'bo', ...
        'MarkerFaceColor','b', 'HandleVisibility','off');
    
    for i = 1:length(w_cross)
        text(0.3, w_cross(i), ...
            sprintf('ω=%.2f\nK=%.2f', w_cross(i), K_cross(i)), ...
            'Color','b');
    end
end

% ====================== 6. 出射角（复极点）+ 入射角（复零点） ======================
% --- 复极点出射角 ---
complex_p = p(abs(imag(p)) > 1e-6);
if ~isempty(complex_p)
    p_upper = complex_p(imag(complex_p) > 0);  % 仅算上半平面，下半共轭对称
    for pk = p_upper
        % 所有零点到该极点的相角和
        sum_z = sum(angle(pk - z));
        % 其他所有极点到该极点的相角和
        sum_p = sum(angle(pk - p(abs(p-pk) > 1e-6)));
        % 出射角归一化
        theta_d = 180 + rad2deg(sum_z - sum_p);
        theta_d = mod(theta_d + 180, 360) - 180;
        
        text(real(pk)+0.8, imag(pk)+0.5, ...
            sprintf('出射角: %.1f°', theta_d), ...
            'Color', [0.6 0 0.6], 'FontSize',9);
    end
end

% --- 复零点入射角 ---
complex_z = z(abs(imag(z)) > 1e-6);
if ~isempty(complex_z)
    z_upper = complex_z(imag(complex_z) > 0);
    for zk = z_upper
        sum_p = sum(angle(zk - p));
        sum_z = sum(angle(zk - z(abs(z-zk) > 1e-6)));
        theta_a = 180 + rad2deg(sum_p - sum_z);
        theta_a = mod(theta_a + 180, 360) - 180;
        
        text(real(zk)+0.8, imag(zk)+0.5, ...
            sprintf('入射角: %.1f°', theta_a), ...
            'Color', [0 0.6 0], 'FontSize',9);
    end
end

legend('Location','best');
hold off;
end

%%
% rlocus_annotated(gs);

%% 比例微分控制环节
% Td = 1;
% cs = 1 + Td*s;
% 最佳阻尼比 
% K = 31.4, Td = 0.08
% 保持 K=20.3
% K = 20.3, Td = 0.1
% 经验估计(阻尼太高，但调节快，无超调)
% K = 50, Td = 0.2
% 只要增大K， 就可以简单地提高响应速度

%% 比例系数的等效根轨迹图
K = 31.4;
ts = K*gs*s/(1+K*gs);
rlocus(ts);

%% 开环增益的根轨迹图
Td = 0.2z;
ks = (1+Td*s)*gs;
rlocus(ks);

%% 时域图像的绘制并标注超调和 2% 调节时间（PD 控制）
% 闭环系统（不含微分项）
closed_ori = feedback(K*gs, 1);
% 闭环系统（含比例-微分项）
closed_pd = feedback(K*(1+Td*s)*gs, 1);

tfinal = 5;        % 仿真总时间（可调整）
t = linspace(0, tfinal, 1000);

% 计算阶跃响应
[y_ori, t_ori] = step(closed_ori, t);
[y_pd,  t_pd]  = step(closed_pd,  t);

% 获取性能指标（超调和 2% 调节时间）
info_pd = stepinfo(y_pd, t_pd, 'SettlingTimeThreshold', 0.02);
os_pd = info_pd.Overshoot;        % 超调百分比
ts_pd = info_pd.SettlingTime;     % 2% 调节时间
peak_val = info_pd.Peak;          % 峰值
peak_time = info_pd.PeakTime;     % 峰值时间

% 判断超调
os_threshold = 1e-2; % 小于此值认为近似0%
has_peak = ~(isempty(peak_time) || peak_time == 0 || os_pd <= os_threshold);

% 绘制在同一图上
figure;
h1 = plot(t_ori, y_ori, 'b-', 'LineWidth', 1.2);
hold on;
h2 = plot(t_pd,  y_pd,  'r-', 'LineWidth', 1.5);

if has_peak
    % 标出峰值点
    [~, idx_peak] = min(abs(t_pd - peak_time));
    xp = t_pd(idx_peak);
    yp = peak_val;
    plot(xp, yp, 'o', 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerSize', 7);
    % 文本偏移
    dx = 0.08 * (t(end) - t(1));      % 水平偏移为时间跨度的 8%
    dy = 0.06 * (max(y_pd) - min(y_pd)); % 垂直偏移为幅值范围的 6%

    txt_x = xp + dx;
    txt_y = yp + dy;
    % 连接线
    plot([xp, txt_x], [yp, txt_y], ':', 'Color', [1 0 0], 'LineWidth', 0.8);
    % 注释文本
    txt1 = sprintf('超调量 = %.2f%%\n峰值时间 = %.3fs', os_pd, peak_time);
    text(txt_x, txt_y, txt1, 'Color', [1 0 0], 'FontSize', 10, 'Interpreter', 'none', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'Margin', 2);
else
    % 无显著超调
    txt_x = t(1) + 0.75*(t(end)-t(1));
    txt_y = min(y_pd) + 0.90*(max(y_pd)-min(y_pd));
    txt1 = '超调量 \approx 0%';
    text(txt_x, txt_y, txt1, 'Color', [1 0 0], 'FontSize', 10, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', 'w', 'Margin', 2);
end

% 调节时间对应点
if ~isempty(ts_pd) && ts_pd > 0 && ts_pd < t(end)
    [~, idx_ts] = min(abs(t_pd - ts_pd));
    xt = t_pd(idx_ts);
    yt = y_pd(idx_ts);
    plot(xt, yt, 's', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 7);

    dx_ts = 0.06 * (t(end) - t(1));
    dy_ts = 0.04 * (max(y_pd) - min(y_pd));
    txt_x2 = xt + dx_ts;
    txt_y2 = yt + dy_ts;
    plot([xt, txt_x2], [yt, txt_y2], ':', 'Color', [0 0 0], 'LineWidth', 0.8);
    txt2 = sprintf('调节时间 (2%%) = %.3fs', ts_pd);
    text(txt_x2, txt_y2, txt2, 'Color', [0 0 0], 'FontSize', 10, 'Interpreter', 'none', ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'Margin', 2);
else
    txt2 = '调节时间 (2%) = N/A';
    text(t(1) + 0.5*(t(end)-t(1)), min(y_pd) + 0.05*(max(y_pd)-min(y_pd)), txt2, ...
         'Color', [0 0 0], 'FontSize', 10, 'Interpreter', 'none', 'HorizontalAlignment', 'center', ...
         'BackgroundColor', 'w', 'Margin', 2);
end
legend([h2,h1], {"增加比例-微分控制后","原系统"}, 'Location', 'Best');

hold off;


