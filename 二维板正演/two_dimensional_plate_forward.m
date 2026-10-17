%% two_dimensional_plate_forward.m
% 二维倾斜板重力与磁力正演
% 对应“二维板正演实验”PPT 的三幅图：模型、重力异常、磁异常。
% 坐标约定：x 向右为正，z 向下为正，观测面 z = 0。
% 几何模型为：上、下边界水平，长度 2L 的两条侧边与 x 轴夹角为 alpha。

clear; clc; close all;

%% 1. 参数（可直接修改，比较模型参数变化）
x0 = 1000.0;                % 板中心 x 坐标 (m)
z0 = 1000.0;                % 板中心 z 坐标 (m)
b  = 200.0;                 % 半宽，PPT: 2b = 400 m
L  = 400.0;                 % 半长，PPT: 2L = 800 m
alpha_deg = 60.0;           % 板倾角 alpha (deg)，可改为 45, 60 等进行比较
i_deg = 60.0;   % 改为 45.0、60.0 等进行比较
M         = 2000.0;         % 磁化强度（按 PPT 给定的数值）
sigma     = 2.67e3;         % 剩余密度 (kg/m^3)，PPT: 2.67 g/cm^3
G         = 6.67430e-11;   % 万有引力常数 (m^3 kg^-1 s^-2)
z_obs     = 0.0;            % 观测面深度 (m)
x_obs     = 0:20:2000;     % PPT: 101 个观测点

alpha = deg2rad(alpha_deg);
iang  = deg2rad(i_deg);

%% 2. 构建二维板模型：A-B-D-C 为顺时针四个角点
% 板长方向的单位向量为 [cos(alpha), sin(alpha)]。
% alpha = 90 deg 时，板为直立有限延深板状体。
u = [cos(alpha), sin(alpha)];
top_mid = [x0, z0] - L * u;
bot_mid = [x0, z0] + L * u;
A = top_mid + [-b, 0];
B = top_mid + [ b, 0];
C = bot_mid + [-b, 0];
D = bot_mid + [ b, 0];

% 与 PPT 的 r1-r4、phi1-phi4 完全对应：A、B、C、D。
corners = [A; B; C; D];

%% 3. 计算重力异常与磁异常分量
% 重力采用二维体元积分，适用于任意倾角的平行四边形板；
% 磁异常采用 PPT 中给出的有限延深二维板解析式。
dg_ms2 = gravity_2d_plate(x_obs, z_obs, x0, z0, b, L, alpha, sigma, G, 100, 200);
[dX, dZ] = plate_magnetic(x_obs, z_obs, corners, alpha, iang, M);
dg_mGal = dg_ms2 * 1e5;  % 1 m/s^2 = 10^5 mGal

% 总场异常：将 (DeltaX, DeltaZ) 投影到磁化方向。
dT = dX .* cos(iang) + dZ .* sin(iang);

mu0 = 4*pi*1e-7;

dX_nT = mu0 * dX * 1e9;
dZ_nT = mu0 * dZ * 1e9;
dT_nT = mu0 * dT * 1e9;

%% 4. 图 1：二维板模型
figure('Color', 'w', 'Position', [100, 120, 850, 500]);
hold on; box on; grid on;
patch([A(1), B(1), D(1), C(1)], [A(2), B(2), D(2), C(2)], ...
      [0.18, 0.45, 0.72], 'FaceAlpha', 0.55, 'EdgeColor', [0.05, 0.15, 0.30], 'LineWidth', 1.8);
plot(x_obs, z_obs + zeros(size(x_obs)), 'k-', 'LineWidth', 1.5);
plot(x0, z0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
text(x0 + 25, z0 + 35, sprintf('Center (%.0f, %.0f) m', x0, z0));
plot(corners(:,1), corners(:,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
text(A(1)-45, A(2)-25, 'A'); text(B(1)+15, B(2)-25, 'B');
text(C(1)-45, C(2)+45, 'C'); text(D(1)+15, D(2)+45, 'D');
xlabel('x (m)'); ylabel('z (m)');
title(sprintf('二维板模型  |  2b = %.0f m, 2L = %.0f m, \\alpha = %.0f^\\circ', ...
    2*b, 2*L, alpha_deg));
axis equal; xlim([min(x_obs), max(x_obs)]); ylim([-100, max([C(2), D(2)]) + 200]);
set(gca, 'YDir', 'reverse');

%% 5. 图 2：重力异常
figure('Color', 'w', 'Position', [130, 150, 850, 430]);
plot(x_obs, dg_mGal, 'b-', 'LineWidth', 2.0); grid on; box on;
xlabel('x (m)'); ylabel('\Delta g (mGal)');
title('二维板重力正演异常');

%% 6. 图 3：磁力异常
figure('Color', 'w', 'Position', [160, 180, 850, 430]);
plot(x_obs, dX, 'Color', [0.85, 0.20, 0.15], 'LineWidth', 1.8); hold on;
plot(x_obs, dZ, 'Color', [0.10, 0.45, 0.75], 'LineWidth', 1.8);
plot(x_obs, dT, 'k-', 'LineWidth', 2.0); grid on; box on;
xlabel('x (m)'); ylabel('磁异常 (nT)');
title(sprintf('二维板磁力正演异常  |  i = %.0f^\\circ, M = %.0f', i_deg, M));
legend('\Delta X', '\Delta Z', '\Delta T', 'Location', 'best');

%% 7. 参数对比示例（取消下行注释即可使用）
% alpha_deg = 45;  i_deg = 45;  % 重新运行脚本，比较倾角和磁偏角的影响。

function [dX, dZ] = plate_magnetic(xk, zk, corners, alpha, iang, M)
% 依据 PPT 中有限延深二维板的解析表达式计算磁异常。
% corners 的行序必须为 [A; B; C; D]；角 phi 从 +x 轴顺时针（即 +z）计算。
    A = corners(1,:); B = corners(2,:); C = corners(3,:); D = corners(4,:);
    n = numel(xk);
    r = zeros(4, n);
    phi = zeros(4, n);
    pts = [A; B; C; D];
    for k = 1:4
        dx = pts(k,1) - xk;
        dz = pts(k,2) - zk;
        r(k,:) = hypot(dx, dz);
        phi(k,:) = atan2(dz, dx);
    end

    % 连续展开相位，避免观测点跨越 atan2 的 pi/-pi 分支时出现伪跳变。
    phi = unwrap(phi, [], 2);
    q = phi(1,:) - phi(2,:) - phi(3,:) + phi(4,:);
    log_ratio = log((r(2,:) .* r(3,:)) ./ (r(1,:) .* r(4,:)));

    dX = M/(2*pi) * sin(alpha) .* ( ...
        log_ratio .* cos(alpha - iang) - sin(alpha - iang) .* q );
    dZ = M/(2*pi) * sin(alpha) .* ( ...
        sin(alpha - iang) .* log_ratio + cos(alpha - iang) .* q );
end

function dg = gravity_2d_plate(xobs, zobs, x0, z0, b, L, alpha, sigma, G, nb, nL)
% 用中点体元积分计算二维无限延伸板的竖直重力异常。
% 单位走向长度上的面积元 dA 的贡献为 2*G*sigma*(z-zobs)/r^2*dA。
% nb、nL 分别为板宽和板长方向剖分数；增加它们可提高精度。
    sb = ((1:nb) - 0.5) / nb * 2*b - b;
    sL = ((1:nL) - 0.5) / nL * 2*L - L;
    [SB, SL] = meshgrid(sb, sL);
    xc = x0 + SL * cos(alpha) + SB;
    zc = z0 + SL * sin(alpha);
    dA = (2*b/nb) * (2*L/nL);

    dg = zeros(size(xobs));
    for k = 1:numel(xobs)
        dx = xc - xobs(k);
        dz = zc - zobs;
        kernel = dz ./ (dx.^2 + dz.^2);
        dg(k) = 2 * G * sigma * sum(kernel(:)) * dA;
    end
end
