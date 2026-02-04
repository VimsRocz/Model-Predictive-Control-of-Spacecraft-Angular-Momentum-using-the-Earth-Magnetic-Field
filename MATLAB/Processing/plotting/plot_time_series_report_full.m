function fig = plot_time_series_report_full(P, S, label)
%PLOT_TIME_SERIES_REPORT_FULL End-to-end report for full attitude+wheels sim.
%
% Required fields in S:
%   t, h_w, m, B_body, tau_ext_body, tau_mtq_body, tau_rw_body, q_ib, q_ib_des, w_body, w_des_body

if nargin < 3
    label = "full";
end

t_s = S.t(:);
h = ensure_nx3(S.h_w);
m = ensure_nx3(S.m);
B_T = ensure_nx3(S.B_body);
tau_ext = ensure_nx3(S.tau_ext_body);
tau_mtq = ensure_nx3(S.tau_mtq_body);
tau_rw = ensure_nx3(S.tau_rw_body);

q = S.q_ib;
qdes = S.q_ib_des;
w = ensure_nx3(S.w_body);
wdes = ensure_nx3(S.w_des_body);

% Derived signals
h_norm = vecnorm(h, 2, 2);
m_norm = vecnorm(m, 2, 2);
B_uT = 1e6 * B_T;
B_norm_uT = vecnorm(B_uT, 2, 2);

tau_des = (-P.KH * h.').';
Bnorm = vecnorm(B_T, 2, 2);
b_hat = B_T ./ max(Bnorm, 1e-30);
tau_par = dot(tau_des, b_hat, 2) .* b_hat;
tau_perp = tau_des - tau_par;

tau_mtq_uNm = 1e6 * vecnorm(tau_mtq,2,2);
tau_ext_uNm = 1e6 * vecnorm(tau_ext,2,2);
tau_rw_uNm = 1e6 * vecnorm(tau_rw,2,2);

mtq_dotB = dot(tau_mtq, b_hat, 2);          % should be ~0
mtq_dotB_uNm = 1e6 * abs(mtq_dotB);

mtq_cos = dot(tau_mtq, b_hat, 2) ./ max(vecnorm(tau_mtq,2,2), 1e-30);
mtq_ang_deg = acosd(max(-1, min(1, mtq_cos))); % angle between tau_mtq and B

% Dipole saturation markers
m_max = P.m_max(:).';
sat_m = (abs(m(:,1)) >= (m_max(1) - 1e-12)) | (abs(m(:,2)) >= (m_max(2) - 1e-12)) | (abs(m(:,3)) >= (m_max(3) - 1e-12));
sat_pct = 100 * mean(sat_m);

% Attitude error angle (deg)
att_err_deg = zeros(size(t_s));
for k = 1:numel(t_s)
    qe = quatmultiply(quatconj(qdes(k,:)), q(k,:));
    if qe(1) < 0; qe = -qe; end
    att_err_deg(k) = 2 * acosd(max(-1, min(1, qe(1))));
end

w_err = w - wdes;
w_err_dps = rad2deg(w_err); %#ok<NASGU>

h0 = h_norm(1);
hf = h_norm(end);
h_red_pct = 100 * (1 - hf / max(h0, 1e-30));

fig = figure('Color','w','Name',sprintf('MTQ Full Report (%s)', label));
fig.Position(3:4) = [1550 1100];
tl = tiledlayout(fig, 6, 2, 'TileSpacing','compact', 'Padding','compact');

title(tl, sprintf('End-to-End MTQ Momentum Desaturation + Attitude (full model) (%s)', label), 'FontWeight','bold');
subtitle(tl, sprintf(['Ts=%.3g s, Tend=%.0f s, bfield=%s | m_{max}=[%.3g %.3g %.3g] A·m^2 | ' ...
    '||h_w||: %.3e → %.3e (%.1f%% reduction) | MTQ sat: %.1f%%'], ...
    P.Ts, P.Tend, bfield_label(P), P.m_max(1), P.m_max(2), P.m_max(3), h0, hf, h_red_pct, sat_pct));

% Wheel momentum components
nexttile(tl);
plot(t_s, h, 'LineWidth', 1.1); grid on;
xlabel('time [s]'); ylabel('h_w [N·m·s]');
legend('h_x','h_y','h_z','Location','best');
title('Reaction Wheel Momentum Components');

% ||h_w||
nexttile(tl);
plot(t_s, h_norm, 'LineWidth', 1.3); grid on;
xlabel('time [s]'); ylabel('||h_w|| [N·m·s]');
title('Wheel Momentum Magnitude');

% Dipole command components
nexttile(tl);
axm = gca;
co = axm.ColorOrder;
hold on; grid on;
for i = 1:3
    plot(t_s, m(:,i), '-', 'Color', co(i,:), 'LineWidth', 1.2);
end
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
legend('m_x','m_y','m_z','Location','best');
title('MTQ Dipole Command (body frame)');
ylim(axm, 1.15*[-max(m_max) max(m_max)]);

% ||m|| + saturation markers
nexttile(tl);
plot(t_s, m_norm, 'LineWidth', 1.3); hold on; grid on;
if any(sat_m)
    scatter(t_s(sat_m), m_norm(sat_m), 12, 'filled', 'MarkerFaceAlpha', 0.5);
    legend('||m||','saturated','Location','best');
else
    legend('||m||','Location','best');
end
xlabel('time [s]'); ylabel('||m|| [A·m^2]');
title('Dipole Magnitude');

% B components (body)
nexttile(tl);
plot(t_s, B_uT, 'LineWidth', 1.1); grid on;
xlabel('time [s]'); ylabel('B [µT]');
legend('B_x','B_y','B_z','Location','best');
title('Earth Magnetic Field in Body Frame');

% ||B||
nexttile(tl);
plot(t_s, B_norm_uT, 'LineWidth', 1.3); grid on;
xlabel('time [s]'); ylabel('||B|| [µT]');
title('Magnetic Field Magnitude');

% Torque magnitudes
nexttile(tl);
plot(t_s, [tau_ext_uNm tau_mtq_uNm tau_rw_uNm], 'LineWidth', 1.1); grid on;
xlabel('time [s]'); ylabel('||\\tau|| [µN·m]');
legend('||\\tau_{ext}||','||\\tau_{mtq}||','||\\tau_{rw}||','Location','best');
title('Torque Magnitudes (body)');

% Perpendicular-to-B limitation
nexttile(tl);
tau_des_uNm = 1e6 * vecnorm(tau_des,2,2);
tau_perp_uNm = 1e6 * vecnorm(tau_perp,2,2);
tau_par_uNm = 1e6 * vecnorm(tau_par,2,2);
plot(t_s, [tau_des_uNm tau_perp_uNm tau_mtq_uNm tau_par_uNm mtq_dotB_uNm], 'LineWidth', 1.1); grid on;
xlabel('time [s]'); ylabel('[µN·m]');
legend('||\\tau_{des}||','||\\tau_{des,\\perp B}||','||\\tau_{mtq}||','||\\tau_{des,\\parallel B}||','|\\tau_{mtq}\\cdot\\hat{b}|','Location','best');
title('MTQ Geometry: \\tau_{mtq} is always \\perp B');

% MTQ angle to B (should be ~90 deg)
nexttile(tl);
plot(t_s, mtq_ang_deg, 'LineWidth', 1.2); grid on;
xlabel('time [s]'); ylabel('angle [deg]');
title('Angle(\\tau_{mtq}, B) (ideal = 90°)');
ylim([0 180]);

% Attitude tracking error
nexttile(tl);
plot(t_s, att_err_deg, 'LineWidth', 1.2); grid on;
xlabel('time [s]'); ylabel('angle error [deg]');
title('Attitude Error to LVLH Reference');

% Equation summary tile
nexttile(tl, [1 2]);
axis off;
set(gca,'XLim',[0 1], 'YLim',[0 1]);
eq = sprintf([ ...
    'Full dynamics (body):\n' ...
    '  q_dot = 0.5 * q ⊗ [0, ω]\n' ...
    '  I*ω_dot = τ_ext + τ_mtq + τ_rw − ω×(Iω + h_w)\n' ...
    '  h_w_dot = −τ_rw\n' ...
    '  τ_mtq = m × B_body   (always ⟂B)\n' ...
    '\nDesaturation (baseline):\n' ...
    '  τ_des = −K_H*h_w\n' ...
    '  m_req = (B × τ_des) / ||B||^2 ,  m_cmd = sat(m_req, ±m_max)\n' ...
    '\nAttitude reference:\n' ...
    '  LVLH (nadir) tracking with wheel PD torque τ_rw.\n' ...
    '\nVisualization hint:\n' ...
    '  The unachievable component is τ_des,∥B; MTQs can only generate τ in the plane ⟂B.\n' ...
    ]);
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

function s = bfield_label(P)
if isfield(P,'env') && isfield(P.env,'bfield_model')
    s = string(P.env.bfield_model);
else
    s = "dipole";
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

