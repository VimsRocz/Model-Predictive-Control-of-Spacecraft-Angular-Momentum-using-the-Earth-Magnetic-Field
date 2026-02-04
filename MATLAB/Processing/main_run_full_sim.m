function main_run_full_sim(desatMode, P)
%MAIN_RUN_FULL_SIM Run full attitude + wheel + MTQ simulation and plot.

clc;

if nargin < 1 || strlength(desatMode) == 0
    desatMode = "baseline";
end

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
matlab_proc = fullfile(repo_root, 'MATLAB', 'Processing');
addpath(matlab_proc);
addpath(fullfile(repo_root, 'MATLAB', 'Input', 'params'));
addpath(fullfile(matlab_proc, 'models'));
addpath(fullfile(matlab_proc, 'controllers'));
addpath(fullfile(matlab_proc, 'plotting'));

if nargin < 2 || isempty(P)
    P = params_default();
end

S = simulate_mtq_full(P, desatMode);

out_dir = fullfile(repo_root, "MATLAB", "Output", "full");
if ~exist(out_dir, "dir")
    mkdir(out_dir);
end

safe_mode = string(regexprep(char(S.desatMode), '[^A-Za-z0-9]+', '_'));
out_file_mode = fullfile(out_dir, "sim_out_full_" + safe_mode + ".mat");
out_file_latest = fullfile(out_dir, "sim_out_full.mat");

save(out_file_mode, "-struct", "S");
save(out_file_mode, "P", "-append"); % keep P at top-level too
try
    copyfile(out_file_mode, out_file_latest, "f");
catch
end

fprintf('MATLAB full output: %s\n', out_file_mode);
fprintf('  ||h_w(0)|| = %.4e [N·m·s]\n', norm(S.h_w(1,:)));
fprintf('  ||h_w(end)|| = %.4e [N·m·s]\n', norm(S.h_w(end,:)));

plot_results_full(out_file_mode);
end
