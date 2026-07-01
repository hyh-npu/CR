%% 原系统
s = tf('s');
K = 10;
g = 4*K/(s*(s+2));

%% 原系统时域特性
closed = feedback(g,1);
[res_ori, t] = step(closed,20);
plot(t, res_ori);

%%
function info = plotFreqCharacteristics(sys, w, figTitle)
% PLOTFREQUENCYCHARACTERISTICS  绘制系统的幅频/相频特性并标注增益/相角裕度
%   info = plotFreqCharacteristics(sys)
%   info = plotFreqCharacteristics(sys, w)
%   info = plotFreqCharacteristics(sys, w, figTitle)
%
% 输入：
%   sys      - LTI 传递函数或模型（tf, zpk, ss）
%   w        - （可选）频率向量（rad/s），默认 logspace(-1,3,1000)
%   figTitle - （可选）图像主标题字符串（显示在幅频子图上方）
%
% 输出 info 结构体包含：
%   info.w        - 使用的频率向量
%   info.mag      - 线性幅值
%   info.phase    - 相位（度）
%   info.mag_dB   - 幅值 (dB)
%   info.GM       - 幅值裕度 (dB)
%   info.PM       - 相角裕度 (deg)
%   info.Wcg      - 增益交叉频率 (rad/s)
%   info.Wcp      - 相角交叉频率 (rad/s)
%   info.Fig      - 图形句柄

if nargin < 2 || isempty(w)
    w = logspace(-1, 3, 1000); % 默认 0.1 ~ 1000 rad/s
end
if nargin < 3 || isempty(figTitle)
    figTitle = '频域特性';
end

% 计算频率响应
[mag, phase] = bode(sys, w);
mag = squeeze(mag);
phase = squeeze(phase);
mag_dB = 20*log10(mag);

% 计算裕度信息（margin 返回以线性比例的增益裕度和频率）
[GM_lin, PM, Wcg, Wcp] = margin(sys);
% 转换 GM 为 dB（若 GM_lin 为 inf，则设为 Inf）
if isfinite(GM_lin) && GM_lin > 0
    GM_dB = 20*log10(GM_lin);
else
    GM_dB = Inf;
end

% 绘图
fig = figure('Color','w');
% 幅频
subplot(2,1,1);
h1 = semilogx(w, mag_dB, 'b-', 'LineWidth', 1.5);
grid on;
ylabel('幅值 (dB)');
title(figTitle);    % 使用自定义标题
hold on;

% 若存在增益交叉频率，标注交叉点和幅值裕度
if ~isempty(Wcg) && isfinite(Wcg) && Wcg>0 && Wcg < max(w)
    % 找到最近的频率索引用于绘图定位
    [~, idx] = min(abs(w - Wcg));
    xg = w(idx);
    yg = mag_dB(idx);
    plot(xg, yg, 'o', 'Color', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 7);
    % 在幅值曲线下方显示 GM（dB）
    txt_x = xg * 10^(0.08); % 右移若干倍
    txt_y = yg - 0.08*(max(mag_dB)-min(mag_dB));
    plot([xg, txt_x], [yg, txt_y], ':', 'Color', [0 0.4470 0.7410], 'LineWidth', 0.8);
    if isfinite(GM_dB)
        txt = sprintf('GM = %.2f dB \n w_{cg} = %.3g rad/s', GM_dB, Wcg);
    else
        txt = sprintf('GM = Inf \n w_{cg} = %.3g rad/s', Wcg);
    end
    text(txt_x, txt_y, txt, 'Color', [0 0.4470 0.7410], 'FontSize', 9, ...
        'HorizontalAlignment','left','VerticalAlignment','bottom', ...
        'BackgroundColor','w','Margin',2);
end
hold off;

% 相频
subplot(2,1,2);
h2 = semilogx(w, phase, 'Color', '#d95319', 'LineWidth', 1.5);
grid on;
xlabel('角频率 \omega (rad/s)');
ylabel('相位 (°)');
hold on;

% 若存在相角交叉频率，标注交叉点和相角裕度
if ~isempty(Wcp) && isfinite(Wcp) && Wcp>0 && Wcp < max(w)
    [~, idx2] = min(abs(w - Wcp));
    xp = w(idx2);
    yp = phase(idx2);
    plot(xp, yp, 's', 'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerSize', 7);
    txt_x2 = xp * 10^(0.06);
    txt_y2 = yp + 0.06*(max(phase)-min(phase));
    plot([xp, txt_x2], [yp, txt_y2], ':', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 0.8);
    txt2 = sprintf('PM = %.2f° \n w_{cp} = %.3g rad/s', PM, Wcp);
    text(txt_x2, txt_y2, txt2, 'Color', [0.8500 0.3250 0.0980], 'FontSize', 9, ...
        'HorizontalAlignment','left','VerticalAlignment','bottom', ...
        'BackgroundColor','w','Margin',2);
end
hold off;

% 返回结构体
info.w = w;
info.mag = mag;
info.phase = phase;
info.mag_dB = mag_dB;
info.GM = GM_dB;
info.PM = PM;
info.Wcg = Wcg;
info.Wcp = Wcp;
info.Fig = fig;
end


%% 原系统频域特性
info = plotFreqCharacteristics(g,[],"原系统开环特性");



%% 串联超前校正
% a = 3.543, T = 0.0861, wc0 = 3.278, wd0 = 11.6144
wc = 3.278; wd = 13.6144; % wc < wd
forward = (1+s/wc)/(1+s/wd);
info  = plotFreqCharacteristics(forward,[],"超前校正装置特性");
%
info = plotFreqCharacteristics(forward*g,[],"超前校正后系统开环特性");

% 比较：绘制校正前后闭环阶跃响应（同一图）
figure('Color','w');
% 校正后（超前）闭环响应
closed_for = feedback(forward*g, 1);
[res_for, t_for] = step(closed_for, 20);
h1 = plot(t_for, res_for, 'b-', 'LineWidth', 1.5);
hold on;
% 原系统闭环响应（之前计算的 res_ori, t）
h2 = plot(t, res_ori, 'r--', 'LineWidth', 1.5);
hold off;
grid on;
xlabel('时间 (s)');
ylabel('输出');
title('超前校正前后闭环阶跃响应比较');
legend([h2, h1], {'原系统', '超前校正后'}, 'Location', 'Best');



%% 迟后-超前校正
% a = 3.3978, w = 4.46323, wx0 = wx, phi = 33.04, 10lga = 5.312
wa = 0.131356; wb = 0.446323; wc = 2.4213; wd = 8.2271; % wa < wb < wc < wd
bcfw = (1+s/wb)*(1+s/wc)/((1+s/wa)*(1+s/wd));
info  = plotFreqCharacteristics(bcfw,[],"迟后超前校正装置特性");
%
info  = plotFreqCharacteristics(bcfw*g,[],"迟后-超前校正后系统开环特性");

% 新增：比较迟后-超前校正前后闭环阶跃响应（同一图）
figure('Color','w');
closed_bcfw = feedback(bcfw*g, 1);
[res_bcfw, t_bcfw] = step(closed_bcfw, 20);
h3 = plot(t_bcfw, res_bcfw, 'b-', 'LineWidth', 1.5);
hold on;
% 原系统闭环响应（之前计算的 res_ori, t）
h4 = plot(t, res_ori, 'r--', 'LineWidth', 1.5);
hold off;
grid on;
xlabel('时间 (s)');
ylabel('输出');
title('迟后-超前校正前后闭环阶跃响应比较');
legend([h4, h3], {'原系统', '迟后-超前校正后'}, 'Location', 'Best');


%% 串联迟后校正
% b = 1041; T = 42.1877; wa0 = 0.0237036; wb0 = 0.2277
wa = 0.0207036; wb = 0.2277; % wa < wb
lag = (1 + s/wb) / (1 + s/wa);
info = plotFreqCharacteristics(lag, [], "串联迟后校正装置特性");
%
info = plotFreqCharacteristics(lag*g, [], "串联迟后校正后系统开环特性");


figure('Color','w');
closed_lag = feedback(lag*g, 1);
[res_lag, t_lag] = step(closed_lag, 20); % 迟后校正通常较慢，延长仿真时间
h5 = plot(t_lag, res_lag, 'b-', 'LineWidth', 1.5);
hold on;

h6 = plot(t, res_ori, 'r--', 'LineWidth', 1.5);
hold off;
grid on;
xlabel('时间 (s)');
ylabel('输出');
title('串联迟后校正前后闭环阶跃响应比较');
legend([h6, h5], {'原系统', '串联迟后校正后'}, 'Location', 'Best');

%% 比较三种校正方式的开环波特图（同一图）
figure('Color','w');
hold on;
w = logspace(-2, 2, 1000);
[mag_ori, ph_ori] = bode(g, w); mag_ori = squeeze(mag_ori); ph_ori = squeeze(ph_ori);
[mag_for, ph_for] = bode(forward*g, w); mag_for = squeeze(mag_for); ph_for = squeeze(ph_for);
[mag_bcfw, ph_bcfw] = bode(bcfw*g, w); mag_bcfw = squeeze(mag_bcfw); ph_bcfw = squeeze(ph_bcfw);
[mag_lag, ph_lag] = bode(lag*g, w); mag_lag = squeeze(mag_lag); ph_lag = squeeze(ph_lag);

subplot(2,1,1);
semilogx(w, 20*log10(mag_ori), 'k--', 'LineWidth', 1);
hold on;
semilogx(w, 20*log10(mag_for), 'b-', 'LineWidth', 1.2);
semilogx(w, 20*log10(mag_bcfw), 'g-', 'LineWidth', 1.2);
semilogx(w, 20*log10(mag_lag), 'm-', 'LineWidth', 1.2);
grid on; ylabel('幅值 (dB)'); title('开环波特图比较'); legend('原系统','串联超前','迟后-超前','串联迟后','Location','Best');

subplot(2,1,2);
semilogx(w, ph_ori, 'k--', 'LineWidth', 1);
hold on;
semilogx(w, ph_for, 'b-', 'LineWidth', 1.2);
semilogx(w, ph_bcfw, 'g-', 'LineWidth', 1.2);
semilogx(w, ph_lag, 'm-', 'LineWidth', 1.2);
grid on; xlabel('\omega (rad/s)'); ylabel('相位 (°)');

hold off;

% 比较三种校正方式的闭环阶跃响应（同一图）
figure('Color','w');
closed_ori = feedback(g,1);
closed_for = feedback(forward*g,1);
closed_bcfw = feedback(bcfw*g,1);
closed_lag = feedback(lag*g,1);

tmax = 10; % 保持原意
t_common = linspace(0, tmax, 1000)'; % 统一时间向量

[y1, ~] = step(closed_ori, t_common);
[y2, ~] = step(closed_for, t_common);
[y3, ~] = step(closed_bcfw, t_common);
[y4, ~] = step(closed_lag, t_common);

% 绘制闭环阶跃响应（时间-输出图）
h_ori = plot(t_common, y1, 'k--', 'LineWidth', 1.2); hold on;
h_for = plot(t_common, y2, 'b-',  'LineWidth', 1.5);
h_bcfw = plot(t_common, y3, 'g-', 'LineWidth', 1.5);
h_lag = plot(t_common, y4, 'm-', 'LineWidth', 1.5);
hold off;
grid on;
xlabel('时间 (s)');
ylabel('输出');
title('校正前后闭环阶跃响应比较（时间-输出图）');
legend([h_ori, h_for, h_bcfw, h_lag], {'原系统','串联超前校正后','迟后-超前校正后','串联迟后校正后'}, 'Location','Best');
