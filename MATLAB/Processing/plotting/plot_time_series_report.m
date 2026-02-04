function fig = plot_time_series_report(P, t_s, x, m, B_T, tau_ext_Nm, tau_mtq_Nm, tau_total_Nm, label)
%PLOT_TIME_SERIES_REPORT End-to-end time-series report with labels + limits.
%
% Inputs are expected as Nx3 arrays (except t_s Nx1).

if nargin < 9
    label = "run";
end

[t_s, x, m, B_T, tau_ext_Nm, tau_mtq_Nm, tau_total_Nm] = normalize_inputs(t_s, x, m, B_T, tau_ext_Nm, tau_mtq_Nm, tau_total_Nm);

% Derived signals
x_norm = vecnorm(x, 2, 2);
m_norm = vecnorm(m, 2, 2);
B_uT = 1e6 * B_T;
B_norm_uT = vecnorm(B_uT, 2, 2);

tau_des = (-P.KH * x.').';
Bnorm2 = sum(B_T.^2, 2);
m_req = cross(B_T, tau_des, 2) ./ max(Bnorm2, 1e-30); % requested (pre-saturation) dipole
% For readability, also compute a clipped version of the requested dipole
% (so the commanded signal isn't visually flattened by the larger request).
m_max = P.m_max(:).';
m_req_clip = min(max(m_req, -m_max), m_max);
Bnorm = vecnorm(B_T, 2, 2);
b_hat = B_T ./ max(Bnorm, 1e-30);
tau_parallel = dot(tau_des, b_hat, 2) .* b_hat;
tau_perp = tau_des - tau_parallel;

% The magnetic torque from MTQ should be perpendicular to B.
tau_mtq_dotB = dot(tau_mtq_Nm, b_hat, 2);

% Saturation markers
sat = (abs(m(:,1)) >= (m_max(1) - 1e-12)) | (abs(m(:,2)) >= (m_max(2) - 1e-12)) | (abs(m(:,3)) >= (m_max(3) - 1e-12));
sat_pct = 100 * mean(sat);

x0 = x_norm(1);
xf = x_norm(end);
x_red_pct = 100 * (1 - xf / max(x0, 1e-30));

fig = figure('Color','w','Name',sprintf('MTQ Desaturation Report (%s)', label));
fig.Position(3:4) = [1500 1050];
tl = tiledlayout(fig, 5, 2, 'TileSpacing','compact', 'Padding','compact');

title(tl, sprintf('End-to-End MTQ Momentum Desaturation (%s)', label), 'FontWeight','bold');
subtitle(tl, sprintf(['Ts=%.3g s, Tend=%.0f s, m_{max}=[%.3g %.3g %.3g] A·m^2, B0=%.1f µT@equator surface | ' ...
    '||x||: %.3e → %.3e (%.1f%% reduction) | sat: %.1f%%'], ...
    P.Ts, P.Tend, P.m_max(1), P.m_max(2), P.m_max(3), 1e6*P.B0_T, x0, xf, x_red_pct, sat_pct));

% x components
nexttile(tl);
plot(t_s, x, 'LineWidth', 1.1); grid on;
xlabel('time [s]'); ylabel('x [N·m·s]');
legend('x_x','x_y','x_z','Location','best');
title('Momentum Error Components');

% ||x||
nexttile(tl);
plot(t_s, x_norm, 'LineWidth', 1.3); grid on;
xlabel('time [s]'); ylabel('||x|| [N·m·s]');
title('Momentum Magnitude');

% m components + limits
nexttile(tl);
axm = gca;
co = axm.ColorOrder;
hold on; grid on;
% Requested (pre-saturation) dipole (clipped, dashed; same axis scale as cmd)
for i = 1:3
    col = 0.65*co(i,:) + 0.35*[1 1 1];
    plot(t_s, m_req_clip(:,i), '--', 'Color', col, 'LineWidth', 0.9, 'HandleVisibility','off');
end
% Commanded (saturated) dipole (solid)
for i = 1:3
    plot(t_s, m(:,i), '-', 'Color', co(i,:), 'LineWidth', 1.2);
end
% Per-axis limits (draw as band if all equal, else individual lines)
if max(abs(diff(m_max))) < 1e-12
    yline(m_max(1), 'k--', 'm_{max}');
    yline(-m_max(1), 'k--', '');
else
    yline(m_max(1), 'k--', 'm_{x,max}'); yline(-m_max(1), 'k--', '');
    yline(m_max(2), 'Color',[0.3 0.3 0.3], 'LineStyle','--', 'Label','m_{y,max}');
    yline(-m_max(2), 'Color',[0.3 0.3 0.3], 'LineStyle','--', 'Label','');
    yline(m_max(3), 'Color',[0.5 0.5 0.5], 'LineStyle','--', 'Label','m_{z,max}');
    yline(-m_max(3), 'Color',[0.5 0.5 0.5], 'LineStyle','--', 'Label','');
end
xlabel('time [s]'); ylabel('m [A·m^2]');
legend('m_x (cmd)','m_y (cmd)','m_z (cmd)','Location','best');
title('Dipole Command (solid) vs Requested (clipped dashed) + saturation limits');
text(0.01, 0.92, 'dashed = requested dipole (clipped to ±m_{max} for readability)', 'Units','normalized', 'FontSize', 9);
ylim(axm, 1.15*[-max(m_max) max(m_max)]);

% ||m|| + saturation markers
nexttile(tl);
plot(t_s, m_norm, 'LineWidth', 1.3); hold on; grid on;
plot(t_s, vecnorm(m_req,2,2), '--', 'Color',[0.25 0.25 0.25], 'LineWidth', 1.1);
if any(sat)
    scatter(t_s(sat), m_norm(sat), 12, 'filled', 'MarkerFaceAlpha', 0.5);
    legend('||m_{cmd}||','||m_{req}||','saturated','Location','best');
else
    legend('||m_{cmd}||','||m_{req}||','Location','best');
end
xlabel('time [s]'); ylabel('||m|| [A·m^2]');
title('Dipole Magnitude (command vs request)');

% B components
nexttile(tl);
plot(t_s, B_uT, 'LineWidth', 1.1); grid on;
xlabel('time [s]'); ylabel('B [µT]');
legend('B_x','B_y','B_z','Location','best');
title('Earth Magnetic Field (body/inertial frame)');

% ||B||
nexttile(tl);
plot(t_s, B_norm_uT, 'LineWidth', 1.3); grid on;
xlabel('time [s]'); ylabel('||B|| [µT]');
title('Magnetic Field Magnitude');

% Torques (ext + MTQ + total) magnitude
nexttile(tl);
tau_ext_uNm = 1e6 * vecnorm(tau_ext_Nm,2,2);
tau_mtq_uNm = 1e6 * vecnorm(tau_mtq_Nm,2,2);
tau_tot_uNm = 1e6 * vecnorm(tau_total_Nm,2,2);
plot(t_s, [tau_ext_uNm tau_mtq_uNm tau_tot_uNm], 'LineWidth', 1.1); grid on;
xlabel('time [s]'); ylabel('||\tau|| [µN·m]');
legend('||\tau_{ext}||','||\tau_{mtq}||','||\tau_{total}||','Location','best');
title('Torque Magnitudes');

% Baseline effect / geometry: unachievable component along B
nexttile(tl);
tau_des_uNm = 1e6 * vecnorm(tau_des,2,2);
tau_perp_uNm = 1e6 * vecnorm(tau_perp,2,2);
tau_par_uNm = 1e6 * vecnorm(tau_parallel,2,2);
mtq_dotB_uNm = 1e6 * abs(tau_mtq_dotB); % should be ~0
plot(t_s, [tau_des_uNm tau_perp_uNm tau_mtq_uNm tau_par_uNm mtq_dotB_uNm], 'LineWidth', 1.1); grid on;
xlabel('time [s]'); ylabel('[µN·m]');
legend('||\tau_{des}||','||\tau_{des,\perp B}||','||\tau_{mtq}||','||\tau_{des,\parallel B}||','|\tau_{mtq}\cdot \hat{b}|','Location','best');
title('Baseline Geometry: MTQ cannot generate torque along B');

% Add a compact equation summary in a dedicated tile (avoid covering plots).
nexttile(tl, [1 2]);
axis off;
set(gca,'XLim',[0 1], 'YLim',[0 1]);
eq = sprintf(['Key equations (momentum-only model):\n' ...
    '  tau_mtq = m × B    (always ⟂B)\n' ...
    '  x_{k+1} = x_k + Ts*(tau_ext + tau_mtq)\n' ...
    '\nBaseline (instantaneous):\n' ...
    '  tau_des = -KH*x\n' ...
    '  m_req   = (B × tau_des) / ||B||^2\n' ...
    '  m_cmd   = sat(m_req, ±m_max)\n' ...
    '\nMPC (uses future B along orbit):\n' ...
    '  minimize Σ (x^T Q x + m^T R m + Δm^T S Δm)\n' ...
    '  subject to x_{i+1} = x_i + Ts*(tau_ext,i - [B_i]_× m_i), |m_i|≤m_max']);
text(0.01, 0.98, eq, 'Interpreter','none', 'VerticalAlignment','top', ...
    'FontName','monospace', 'FontSize', 9);

% Remove axes toolbars (prevents exportgraphics from capturing icons)
axs = findall(fig, 'Type','axes');
for i = 1:numel(axs)
    try
        axs(i).Toolbar.Visible = 'off';
    catch
    end
end
end

function [t_s, x, m, B_T, tau_ext_Nm, tau_mtq_Nm, tau_total_Nm] = normalize_inputs(t_s, x, m, B_T, tau_ext_Nm, tau_mtq_Nm, tau_total_Nm)
% Normalize to Nx3 arrays.

t_s = t_s(:);

x = ensure_nx3(x);
m = ensure_nx3(m);
B_T = ensure_nx3(B_T);
tau_ext_Nm = ensure_nx3(tau_ext_Nm);
tau_mtq_Nm = ensure_nx3(tau_mtq_Nm);
tau_total_Nm = ensure_nx3(tau_total_Nm);

n = numel(t_s);
if size(x,1) ~= n
    error('Time vector length (%d) does not match x rows (%d).', n, size(x,1));
end
end

function A = ensure_nx3(A)
if isempty(A)
    A = zeros(0,3);
    return;
end
if size(A,2) == 3
    return;
end
if size(A,1) == 3
    A = A.';
    return;
end
error('Expected Nx3 or 3xN array. Got %dx%d.', size(A,1), size(A,2));
end
