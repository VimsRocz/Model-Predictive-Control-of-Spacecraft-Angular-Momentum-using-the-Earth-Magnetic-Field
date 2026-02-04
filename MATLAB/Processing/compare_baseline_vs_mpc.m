function compare_baseline_vs_mpc(P)
%COMPARE_BASELINE_VS_MPC Compare MATLAB baseline vs MPC runs and plot metrics.
%
% Expects:
%   MATLAB/Output/momentum/sim_out_matlab_baseline.mat
%   MATLAB/Output/momentum/sim_out_matlab_mpc.mat
%
% If P is provided, it is used only for labeling (files are loaded from disk).

repo_root = fileparts(fileparts(mfilename('fullpath')));
matlab_proc = fullfile(repo_root, 'MATLAB', 'Processing');
addpath(matlab_proc);
addpath(fullfile(repo_root, 'MATLAB', 'Input', 'params'));
addpath(fullfile(matlab_proc, 'models'));
addpath(fullfile(matlab_proc, 'controllers'));
addpath(fullfile(matlab_proc, 'plotting'));
out_dir = fullfile(repo_root, 'MATLAB', 'Output', 'momentum');

f_base = fullfile(out_dir, 'sim_out_matlab_baseline.mat');
f_mpc  = fullfile(out_dir, 'sim_out_matlab_mpc.mat');

if ~exist(f_base,'file')
    error("Missing baseline output: %s", f_base);
end
if ~exist(f_mpc,'file')
    error("Missing MPC output: %s", f_mpc);
end

Sb = load(f_base);
Sm = load(f_mpc);

Pb = Sb.P;
Pm = Sm.P;

% Prefer provided P label; otherwise use baseline P
if nargin < 1 || isempty(P)
    P = Pb;
end

tb = (0:size(Sb.x,2)-1) * Pb.Ts;
tm = (0:size(Sm.x,2)-1) * Pm.Ts;

xb = Sb.x.'; % Nx3
xm = Sm.x.'; % Nx3
mb = Sb.m.'; % Nx3
mm = Sm.m.'; % Nx3

xbn = vecnorm(xb,2,2);
xmn = vecnorm(xm,2,2);
mbn = vecnorm(mb,2,2);
mmn = vecnorm(mm,2,2);

satb = sat_pct(Sb.m, Pb.m_max);
satm = sat_pct(Sm.m, Pm.m_max);

fig = figure('Color','w','Name','Baseline vs MPC (MATLAB)');
fig.Position(3:4) = [1400 800];
tl = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');

title(tl, 'Baseline vs MPC (MATLAB)', 'FontWeight','bold');
subtitle(tl, sprintf(['Orbit alt=%.0f km, inc=%.1f deg | B0=%.1f µT | m_{max}=%.2f A·m^2\n' ...
    'Env: bias=%d, periodic=%d, GG=%d, drag=%d'], ...
    P.orbit.alt_m/1e3, P.orbit.inc_deg, 1e6*P.B0_T, P.m_max(1), ...
    truthy(P.env.enable_bias), truthy(P.env.enable_periodic), truthy(P.env.enable_gravity_gradient), truthy(P.env.enable_drag)));

nexttile(tl);
plot(tb, xbn, 'LineWidth', 1.3); hold on;
plot(tm, xmn, '--', 'LineWidth', 1.5);
grid on;
xlabel('time [s]'); ylabel('||x|| [N·m·s]');
legend('baseline','mpc','Location','best');
title(sprintf('Momentum Magnitude (final: %.3e → %.3e)', xbn(end), xmn(end)));

nexttile(tl);
plot(tb, mbn, 'LineWidth', 1.3); hold on;
plot(tm, mmn, '--', 'LineWidth', 1.5);
grid on;
xlabel('time [s]'); ylabel('||m|| [A·m^2]');
legend('baseline','mpc','Location','best');
title(sprintf('Dipole Magnitude (sat%%: %.1f vs %.1f)', satb, satm));

nexttile(tl);
plot(tb, xb, 'LineWidth', 1.0); grid on;
xlabel('time [s]'); ylabel('x [N·m·s]');
legend('x_x','x_y','x_z','Location','best');
title('Baseline x Components');

nexttile(tl);
plot(tm, xm, 'LineWidth', 1.0); grid on;
xlabel('time [s]'); ylabel('x [N·m·s]');
legend('x_x','x_y','x_z','Location','best');
title('MPC x Components');

% Remove axes toolbars (prevents exportgraphics from capturing icons)
axs = findall(fig, 'Type','axes');
for i = 1:numel(axs)
    try
        axs(i).Toolbar.Visible = 'off';
    catch
    end
end

% Report solve time if present
if isfield(Sm, 'solve_time_s') && any(Sm.solve_time_s > 0)
    st = Sm.solve_time_s(:);
    st_ms = 1e3 * st(st > 0);
    if ~isempty(st_ms)
        fprintf('MPC solve time: mean %.2f ms, p95 %.2f ms, max %.2f ms (N=%d)\n', ...
            mean(st_ms), prctile(st_ms,95), max(st_ms), numel(st_ms));
    end
end

out_png = fullfile(out_dir, 'compare_baseline_vs_mpc.png');
exportgraphics(fig, out_png, 'Resolution', 250);
fprintf('Saved: %s\n', out_png);
end

function pct = sat_pct(m, m_max)
m = m.';
mx = m_max(:).';
sat = (abs(m(:,1)) >= (mx(1) - 1e-12)) | (abs(m(:,2)) >= (mx(2) - 1e-12)) | (abs(m(:,3)) >= (mx(3) - 1e-12));
pct = 100 * mean(sat);
end

function v = truthy(x)
v = 0;
try
    v = double(logical(x));
catch
end
end
