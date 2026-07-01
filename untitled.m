clear; clc; close all;

%% ========== 参数设置 ==========
% 第一组：固定自然频率，改变阻尼比
wn_fixed = 10;          % 固定ωn = 10 rad/s
zeta_set = [0.3, 0.6, 0.9]; % 三组阻尼比
colors = {'#0072BD', '#D95319', '#EDB120'}; % 对应配色

% 第二组：固定阻尼比，改变自然频率
zeta_fixed = 0.6;       % 固定ζ = 0.6
wn_set = [5, 10, 15];   % 三组自然频率

%% ========== 创建画布：2×2 子图 ==========
figure('Color','w','Position',[100,100,1050,720]);
sgtitle('阻尼比 \zeta 与自然频率 \omega_n 的几何意义及时域对应关系','FontSize',14);

%% ---------- 子图1：固定ωn，不同ζ的s平面极点分布 ----------
subplot(2,2,1);
hold on; grid on; axis equal;
xlabel('实轴 Re(s)'); ylabel('虚轴 Im(s)');
title(['固定 \omega_n = ',num2str(wn_fixed),' rad/s，极点位置']);

% 绘制参考圆弧：半径为ωn，直观体现"极点到原点距离=ωn"
theta_arc = linspace(0, pi, 200);
plot(wn_fixed*cos(theta_arc), wn_fixed*sin(theta_arc), ...
    'k--', 'LineWidth',0.8, 'HandleVisibility','off');

for i = 1:length(zeta_set)
    zeta = zeta_set(i);
    sigma = zeta * wn_fixed;        % 极点实部绝对值
    wd = wn_fixed * sqrt(1-zeta^2); % 阻尼自然频率（虚部）
    p_upper = -sigma + 1j*wd;       % 上半平面极点
    
    % 绘制极点（叉号）
    plot(real(p_upper), imag(p_upper), 'x', ...
        'Color',colors{i}, 'LineWidth',1.5, 'MarkerSize',10);
    plot(real(p_upper), -imag(p_upper), 'x', ...
        'Color',colors{i}, 'LineWidth',1.5, 'MarkerSize',10, 'HandleVisibility','off');
    
    % 绘制原点到极点的连线
    plot([0, real(p_upper)], [0, imag(p_upper)], ...
        '-', 'Color',colors{i}, 'LineWidth',1);
    
    % 标注阻尼比
    text(real(p_upper)+0.4, imag(p_upper)+0.4, ...
        ['\zeta = ',num2str(zeta)], 'Color',colors{i}, 'FontSize',9);
end

% 标注几何结论
text(-wn_fixed*0.55, wn_fixed*0.85, ...
    {'极点到原点距离 = \omega_n', 'cos\theta = \zeta'}, ...
    'Color','k', 'FontSize',9);
xlim([-wn_fixed*1.2, wn_fixed*0.2]);
ylim([0, wn_fixed*1.2]);

%% ---------- 子图2：固定ωn，不同ζ的阶跃响应 ----------
subplot(2,2,2);
hold on; grid on;
xlabel('时间 t (s)'); ylabel('输出 c(t)');
title(['固定 \omega_n = ',num2str(wn_fixed),' rad/s，阶跃响应']);
t = linspace(0, 1.5, 500);

for i = 1:length(zeta_set)
    zeta = zeta_set(i);
    num = wn_fixed^2;
    den = [1, 2*zeta*wn_fixed, wn_fixed^2];
    sys = tf(num, den);
    y = step(sys, t);
    plot(t, y, 'Color',colors{i}, 'LineWidth',1.2, ...
        'DisplayName',['\zeta = ',num2str(zeta)]);
end
legend('Location','southeast');
ylim([0, 1.6]);

%% ---------- 子图3：固定ζ，不同ωn的s平面极点分布 ----------
subplot(2,2,3);
hold on; grid on; axis equal;
xlabel('实轴 Re(s)'); ylabel('虚轴 Im(s)');
title(['固定 \zeta = ',num2str(zeta_fixed),'，极点位置']);

% 绘制参考射线：直观体现"同射线ζ相同"
theta_zeta = acos(zeta_fixed);  % 极点与负实轴的夹角
t_ray = linspace(0, max(wn_set)*1.2, 100);
plot(-t_ray*cos(theta_zeta), t_ray*sin(theta_zeta), ...
    'k--', 'LineWidth',0.8, 'HandleVisibility','off');

for i = 1:length(wn_set)
    wn = wn_set(i);
    sigma = zeta_fixed * wn;
    wd = wn * sqrt(1-zeta_fixed^2);
    p_upper = -sigma + 1j*wd;
    
    plot(real(p_upper), imag(p_upper), 'x', ...
        'Color',colors{i}, 'LineWidth',1.5, 'MarkerSize',10);
    plot(real(p_upper), -imag(p_upper), 'x', ...
        'Color',colors{i}, 'LineWidth',1.5, 'MarkerSize',10, 'HandleVisibility','off');
    
    plot([0, real(p_upper)], [0, imag(p_upper)], ...
        '-', 'Color',colors{i}, 'LineWidth',1);
    
    text(real(p_upper)+0.4, imag(p_upper)+0.4, ...
        ['\omega_n = ',num2str(wn)], 'Color',colors{i}, 'FontSize',9);
end

text(-max(wn_set)*0.65, max(wn_set)*0.75, ...
    {'同射线 → \zeta 相同', '到原点距离 = \omega_n'}, ...
    'Color','k', 'FontSize',9);
xlim([-max(wn_set)*1.2, max(wn_set)*0.2]);
ylim([0, max(wn_set)*1.2]);

%% ---------- 子图4：固定ζ，不同ωn的阶跃响应 ----------
subplot(2,2,4);
hold on; grid on;
xlabel('时间 t (s)'); ylabel('输出 c(t)');
title(['固定 \zeta = ',num2str(zeta_fixed),'，阶跃响应']);
t = linspace(0, 2, 500);

for i = 1:length(wn_set)
    wn = wn_set(i);
    num = wn^2;
    den = [1, 2*zeta_fixed*wn, wn^2];
    sys = tf(num, den);
    y = step(sys, t);
    plot(t, y, 'Color',colors{i}, 'LineWidth',1.2, ...
        'DisplayName',['\omega_n = ',num2str(wn)]);
end
legend('Location','southeast');
ylim([0, 1.6]);