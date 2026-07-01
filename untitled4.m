clear; clc;
syms K T positive

%% ========== 1. 符号计算 I_K、I_T ==========
d2 = 24 + K*T;
d1 = 104 + K + 20*K*T;
d0 = 20*K;

n0 = 104; n1 = 24; n2 = 1;

A = [0, 1, 0; 0, 0, 1; -d0, -d1, -d2];
C = [n0, n1, n2];

syms p11 p12 p13 p22 p23 p33
P = [p11, p12, p13; p12, p22, p23; p13, p23, p33];
eq = A.'*P + P*A == -C.'*C;
sol = solve(eq, [p11,p12,p13,p22,p23,p33]);
P_val = [sol.p11, sol.p12, sol.p13; sol.p12, sol.p22, sol.p23; sol.p13, sol.p23, sol.p33];

ISE = simplify(C * P_val * C.');
I_K = simplify(diff(ISE, K) / 2);
I_T = simplify(diff(ISE, T) / 2);

% 转为数值函数句柄
f_IK = matlabFunction(I_K, 'Vars', [K, T]);
f_IT = matlabFunction(I_T, 'Vars', [K, T]);

%% ========== 2. 生成正参数网格，计算数值分布 ==========
T_min = 0.001;   T_max = 0.2;
K_min = 0.1;     K_max = 500;
grid_num = 80;

T_grid = linspace(T_min, T_max, grid_num);
K_grid = linspace(K_min, K_max, grid_num);
[T_mat, K_mat] = meshgrid(T_grid, K_grid);

IK_map = zeros(size(T_mat));
IT_map = zeros(size(T_mat));

for i = 1:numel(T_mat)
    IK_map(i) = f_IK(K_mat(i), T_mat(i));
    IT_map(i) = f_IT(K_mat(i), T_mat(i));
end

%% ========== 3. 同一张图绘制两条隐函数曲线 ==========
figure('Color','w','Position',[100,100,700,550]);

% 绘制 I_K = 0 零值线（红色虚线，正范围内无轨迹）
contour(T_grid, K_grid, IK_map, [0 0], ...
    'r--', 'LineWidth',1.5, 'DisplayName','I_K = 0');
hold on;

% 绘制 I_T = 0 零值线（蓝色实线）
contour(T_grid, K_grid, IT_map, [0 0], ...
    'b-', 'LineWidth',1.6, 'DisplayName','I_T = 0');

xlabel('微分系数 T');
ylabel('开环增益 K');

legend('Location','southeast');
grid on;
axis([T_min, T_max, K_min, K_max]);


%% ========== 4. 输出示例点验证 ==========
% 找K=100对应的最优T
[~, idx_K100] = min(abs(K_grid - 100));
T_opt_100 = fzero(@(t) f_IT(100, t), 0.05);
fprintf('K = 100 时，最优 T = %.4f s\n', T_opt_100);
fprintf('K = 300 时，最优 T = %.4f s\n', fzero(@(t) f_IT(300, t), 0.08));