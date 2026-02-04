function run_simulink_full(P)
%RUN_SIMULINK_FULL Build and run the full-plant Simulink model and plot.

clc;

repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);
addpath(fullfile(repo_root, 'params'));
addpath(fullfile(repo_root, 'models'));
addpath(fullfile(repo_root, 'controllers'));
addpath(fullfile(repo_root, 'plotting'));

close all;
try
    bdclose('all');
catch
end

if nargin < 1 || isempty(P)
    P = params_default();
end
assignin('base','P',P);

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

out_dir = fullfile(repo_root, 'outputs', 'simulink_full');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

out_file = fullfile(out_dir, 'sim_out_full_simulink.mat');
save(out_file, "-struct", "S");
save(out_file, "P", "-append");

fprintf('Simulink full output: %s\n', out_file);

plot_results_full_simulink(out_file);
end

function A = normalize_nx3(A)
% Ensure A is Nx3
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
% Ensure A is Nx4
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
