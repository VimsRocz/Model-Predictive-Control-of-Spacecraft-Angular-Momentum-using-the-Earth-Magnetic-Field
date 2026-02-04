function run_thesis_benchmarks()
%RUN_THESIS_BENCHMARKS Run a small scenario suite and compare baseline vs MPC (MATLAB).
%
% Outputs:
%   outputs/analysis/benchmark_summary.csv
%   outputs/analysis/benchmark_summary.mat
%   outputs/analysis/benchmark_summary.png

repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);
addpath(fullfile(repo_root, 'params'));
addpath(fullfile(repo_root, 'models'));
addpath(fullfile(repo_root, 'controllers'));
addpath(fullfile(repo_root, 'plotting'));

P0 = params_default();

out_dir = fullfile(repo_root, 'outputs', 'analysis');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Scenario definitions (extend as needed)
scenarios = { ...
    struct('name',"Bias+Periodic (default)", 'apply', @(P) P), ...
    struct('name',"Gravity gradient (with bias)", 'apply', @(P) set_fields(P, {'env.enable_gravity_gradient',true, 'env.enable_periodic',false})), ...
    struct('name',"Aerodynamic drag (with bias)", 'apply', @(P) set_fields(P, {'env.enable_drag',true, 'env.enable_periodic',false})), ...
    struct('name',"GG + Drag (with bias)", 'apply', @(P) set_fields(P, {'env.enable_gravity_gradient',true, 'env.enable_drag',true, 'env.enable_periodic',false})), ...
    struct('name',"Gravity gradient only", 'apply', @(P) set_fields(P, {'env.enable_bias',false, 'env.enable_periodic',false, 'env.enable_gravity_gradient',true})), ...
    struct('name',"Aerodynamic drag only", 'apply', @(P) set_fields(P, {'env.enable_bias',false, 'env.enable_periodic',false, 'env.enable_drag',true})), ...
    struct('name',"Low authority (m_max=0.05)", 'apply', @(P) set_fields(P, {'m_max',0.05*ones(3,1)})), ...
    struct('name',"High authority (m_max=0.5)", 'apply', @(P) set_fields(P, {'m_max',0.5*ones(3,1)})), ...
    };

controllers = ["baseline","mpc"];

rows = [];

for si = 1:numel(scenarios)
    sc = scenarios{si};
    P = sc.apply(P0);
    P.N = round(P.Tend / P.Ts);

    for ci = 1:numel(controllers)
        mode = controllers(ci);
        S = simulate_mtq(P, mode);
        M = metrics(S);

        row = struct();
        row.scenario = sc.name;
        row.controller = mode;
        row.x0 = M.x0;
        row.xf = M.xf;
        row.reduction_pct = M.reduction_pct;
        row.int_x = M.int_x;
        row.sat_pct = M.sat_pct;
        row.mpc_mean_ms = M.mpc_mean_ms;
        row.mpc_p95_ms = M.mpc_p95_ms;
        row.mpc_max_ms = M.mpc_max_ms;

        rows = [rows; row]; %#ok<AGROW>
        fprintf('[%s] %s: ||x|| %.3e -> %.3e (%.1f%%) | sat %.1f%%\n', ...
            sc.name, mode, M.x0, M.xf, M.reduction_pct, M.sat_pct);
    end
end

T = struct2table(rows);
csv = fullfile(out_dir, 'benchmark_summary.csv');
mat = fullfile(out_dir, 'benchmark_summary.mat');
writetable(T, csv);
save(mat, 'T', 'scenarios', 'controllers');

% Simple bar chart: final ||x|| per scenario
fig = figure('Color','w','Name','Benchmark Summary');
fig.Position(3:4) = [1400 650];
tl = tiledlayout(fig, 1, 1);
ax = nexttile(tl);

names = unique(T.scenario, 'stable');
xf_base = nan(numel(names),1);
xf_mpc  = nan(numel(names),1);
for i = 1:numel(names)
    r1 = T(strcmp(T.scenario,names{i}) & strcmp(T.controller,'baseline'),:);
    r2 = T(strcmp(T.scenario,names{i}) & strcmp(T.controller,'mpc'),:);
    if ~isempty(r1); xf_base(i) = r1.xf(1); end
    if ~isempty(r2); xf_mpc(i) = r2.xf(1); end
end

bar(ax, [xf_base xf_mpc]);
grid(ax,'on');
ax.XTickLabel = names;
ax.XTickLabelRotation = 20;
ylabel(ax,'final ||x|| [N·m·s]');
legend(ax, {'baseline','mpc'}, 'Location','best');
title(ax,'Final Momentum Magnitude by Scenario (MATLAB)');

png = fullfile(out_dir, 'benchmark_summary.png');
exportgraphics(fig, png, 'Resolution', 250);

fprintf('Wrote: %s\n', csv);
fprintf('Wrote: %s\n', mat);
fprintf('Wrote: %s\n', png);
end

function M = metrics(S)
t = S.t(:);
x = S.x;
m = S.m;

x_norm = vecnorm(x,2,1).';
m_norm = vecnorm(m,2,1).';

M.x0 = x_norm(1);
M.xf = x_norm(end);
M.reduction_pct = 100 * (1 - M.xf / max(M.x0, 1e-30));
M.int_x = trapz(t, x_norm);

m_max = S.P.m_max(:).';
mm = m.';
sat = (abs(mm(:,1)) >= (m_max(1) - 1e-12)) | (abs(mm(:,2)) >= (m_max(2) - 1e-12)) | (abs(mm(:,3)) >= (m_max(3) - 1e-12));
M.sat_pct = 100 * mean(sat);

M.mpc_mean_ms = NaN;
M.mpc_p95_ms = NaN;
M.mpc_max_ms = NaN;
if isfield(S,'solve_time_s') && any(S.solve_time_s > 0)
    st = S.solve_time_s(:);
    st_ms = 1e3 * st(st > 0);
    if ~isempty(st_ms)
        M.mpc_mean_ms = mean(st_ms);
        M.mpc_p95_ms = prctile(st_ms,95);
        M.mpc_max_ms = max(st_ms);
    end
end
end

function P = set_fields(P, kv)
%SET_FIELDS Apply { 'path', value, 'path', value, ... } updates to struct P.
for i = 1:2:numel(kv)
    path = string(kv{i});
    val = kv{i+1};
    parts = split(path,'.');
    switch numel(parts)
        case 1
            P.(parts{1}) = val;
        case 2
            P.(parts{1}).(parts{2}) = val;
        otherwise
            error('Unsupported path: %s', path);
    end
end
end
