function main_run_sim(ctrlMode, P)
clc;

if nargin < 1 || strlength(ctrlMode) == 0
    % Keep MATLAB and Simulink results consistent by default.
    ctrlMode = "baseline";
end

repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);
addpath(fullfile(repo_root, 'MATLAB', 'Input', 'params'));
addpath(fullfile(repo_root, 'models'));
addpath(fullfile(repo_root, 'controllers'));
addpath(fullfile(repo_root, 'plotting'));

if nargin < 2 || isempty(P)
    P = params_default();
end

S = simulate_mtq(P, ctrlMode);

x = S.x;
m = S.m;
tau_ext_log = S.tau_ext_log;
tau_mtq_log = S.tau_mtq_log;
tau_total_log = S.tau_total_log;
B_log = S.B_log;
ctrlMode = S.ctrlMode;
solve_time_s = S.solve_time_s; %#ok<NASGU>

out_dir = fullfile(repo_root, "MATLAB", "Output", "momentum");
if ~exist(out_dir, "dir")
    mkdir(out_dir);
end
safe_mode = string(regexprep(char(ctrlMode), '[^A-Za-z0-9]+', '_'));
out_file_mode = fullfile(out_dir, "sim_out_matlab_" + safe_mode + ".mat");
out_file_latest = fullfile(out_dir, "sim_out_matlab.mat");

save(out_file_mode, "P","x","m","tau_ext_log","tau_mtq_log","tau_total_log","B_log","ctrlMode","solve_time_s");
try
    copyfile(out_file_mode, out_file_latest, "f");
catch
end

fprintf('MATLAB output: %s\n', out_file_mode);
fprintf('  x(0) = [% .4e % .4e % .4e]\n', x(1,1), x(2,1), x(3,1));
fprintf('  m(0) = [% .4e % .4e % .4e]\n', m(1,1), m(2,1), m(3,1));

plot_results(out_file_mode);
end
