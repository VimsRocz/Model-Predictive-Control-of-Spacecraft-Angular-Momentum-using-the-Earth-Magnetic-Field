function run_simulink(P)
%RUN_SIMULINK Build and run the Simulink model, save outputs, and plot.

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

out_dir = fullfile(repo_root, 'outputs', 'simulink');
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
