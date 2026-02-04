function run_simulink(P)
%RUN_SIMULINK Single entry-point for Simulink simulation (momentum-only or full).
%
% Uses:
%   P.plant.model = "momentum" | "full"
%   P.simulink.controller = "baseline" | "mpc"

clc;

repo_root = fileparts(mfilename('fullpath'));
matlab_proc = fullfile(repo_root, 'MATLAB', 'Processing');
simulink_proc = fullfile(repo_root, 'Simulink', 'Processing');
addpath(matlab_proc);
addpath(simulink_proc);
addpath(fullfile(repo_root, 'MATLAB', 'Input', 'params'));
addpath(fullfile(matlab_proc, 'models'));
addpath(fullfile(matlab_proc, 'controllers'));
addpath(fullfile(matlab_proc, 'plotting'));

close all;
try
    bdclose('all');
catch
end

if nargin < 1 || isempty(P)
    P = params_default();
end

if ~isfield(P,'plant') || ~isfield(P.plant,'model')
    P.plant.model = "full";
end
if ~isfield(P,'simulink') || ~isfield(P.simulink,'controller') || strlength(P.simulink.controller) == 0
    if isfield(P,'controller')
        P.simulink.controller = string(P.controller);
    else
        P.simulink.controller = "baseline";
    end
end

% Force Simulink IMF wrappers to reload parameters when the UI changes.
P.simulink.run_id = generate_run_id();

assignin('base','P',P);

plant = lower(string(P.plant.model));

if plant == "full"
    run_simulink_full_internal(P);
else
    run_simulink_momentum_internal(P);
end
end

function run_simulink_momentum_internal(P)
% Momentum-only Simulink run

% If using IGRF, build an Interpreted MATLAB Function version of the model so
% we can call Aerospace Toolbox functions (igrfmagm) without codegen limits.
use_imf = isfield(P,'env') && isfield(P.env,'bfield_model') && strcmpi(string(P.env.bfield_model), "igrf");
if use_imf
    build_simulink_model_imf();
else
    build_simulink_model();
end

simOut = sim('mtq_desat_mpc_model', 'ReturnWorkspaceOutputs','on');

x_log = simOut.get('x_log');
m_log = simOut.get('m_log');
B_log = simOut.get('B_log');
tau_ext_log = simOut.get('tau_ext_log');
tau_mtq_log = simOut.get('tau_mtq_log');
tau_total_log = simOut.get('tau_total_log');

out_dir = fullfile(repo_root, 'Simulink', 'Output', 'momentum');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
out_file = fullfile(out_dir, 'sim_out_simulink.mat');

save(out_file, 'P', 'x_log', 'm_log', 'B_log', 'tau_ext_log', 'tau_mtq_log', 'tau_total_log');

% Quick sanity print to avoid plotting stale/incorrect runs.
x0_sim = x_log.signals.values(1,:).';
m0_sim = m_log.signals.values(1,:).';
fprintf('Simulink output: %s\n', out_file);
fprintf('  x(0) = [% .4e % .4e % .4e]\n', x0_sim(1), x0_sim(2), x0_sim(3));
fprintf('  m(0) = [% .4e % .4e % .4e]\n', m0_sim(1), m0_sim(2), m0_sim(3));

plot_results_simulink(out_file);
end

function run_simulink_full_internal(P)
% Full-plant Simulink run
repo_root = fileparts(mfilename('fullpath'));
matlab_proc = fullfile(repo_root, 'MATLAB', 'Processing');
simulink_proc = fullfile(repo_root, 'Simulink', 'Processing');
addpath(matlab_proc);
addpath(simulink_proc);
addpath(fullfile(repo_root, 'MATLAB', 'Input', 'params'));
addpath(fullfile(matlab_proc, 'models'));
addpath(fullfile(matlab_proc, 'controllers'));
addpath(fullfile(matlab_proc, 'plotting'));

build_simulink_full_model_imf();

simOut = sim('mtq_full_model', 'ReturnWorkspaceOutputs','on');

S = struct();
S.P = P;
S.t = simOut.get('r_eci_log').time(:);

S.r_eci = normalize_nx3(simOut.get('r_eci_log').signals.values);
S.B_eci = normalize_nx3(simOut.get('B_eci_log').signals.values);
S.B_body = normalize_nx3(simOut.get('B_body_log').signals.values);
S.h_w = normalize_nx3(simOut.get('h_w_log').signals.values);
S.m = normalize_nx3(simOut.get('m_log').signals.values);
S.tau_ext_body = normalize_nx3(simOut.get('tau_ext_log').signals.values);
S.tau_mtq_body = normalize_nx3(simOut.get('tau_mtq_log').signals.values);
S.tau_rw_body = normalize_nx3(simOut.get('tau_rw_log').signals.values);
S.q_ib = normalize_nx4(simOut.get('q_ib_log').signals.values);
S.q_ib_des = normalize_nx4(simOut.get('q_ib_des_log').signals.values);
S.w_body = normalize_nx3(simOut.get('w_body_log').signals.values);
S.w_des_body = normalize_nx3(simOut.get('w_des_body_log').signals.values);

S.desatMode = "simulink/" + string(P.simulink.controller);

out_dir = fullfile(repo_root, 'Simulink', 'Output', 'full');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

out_file = fullfile(out_dir, 'sim_out_full_simulink.mat');
save(out_file, "-struct", "S");
save(out_file, "P", "-append");

fprintf('Simulink full output: %s\n', out_file);

plot_results_full_simulink(out_file);

% If a MATLAB run exists, print a quick consistency check.
try
    fmat = fullfile(repo_root, 'MATLAB', 'Output', 'full', 'sim_out_full.mat');
    if exist(fmat,'file')
        Sm = load(fmat);
        t = Sm.t(:);
        hm = Sm.h_w;
        hs = S.h_w;
        ts = S.t(:);
        hs_i = interp1(ts, hs, t, 'linear', 'extrap');
        diff = vecnorm(hm - hs_i, 2, 2);
        fprintf('MATLAB vs Simulink (full) max |h_w diff| = %.3e\n', max(diff));
    end
catch
end
end

function A = normalize_nx3(A)
A = squeeze(A);
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
error('Expected Nx3 or 3xN array after squeeze. Got %dx%d.', size(A,1), size(A,2));
end

function A = normalize_nx4(A)
A = squeeze(A);
if isempty(A)
    A = zeros(0,4);
    return;
end
if size(A,2) == 4
    return;
end
if size(A,1) == 4
    A = A.';
    return;
end
error('Expected Nx4 or 4xN array after squeeze. Got %dx%d.', size(A,1), size(A,2));
end

function id = generate_run_id()
try
    id = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
catch
    id = char(num2str(now, 16));
end
end
