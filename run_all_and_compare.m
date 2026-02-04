function run_all_and_compare(P)
%RUN_ALL_AND_COMPARE Re-run both MATLAB and Simulink flows and compare outputs.

repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);
addpath(fullfile(repo_root, 'MATLAB', 'Input', 'params'));
addpath(fullfile(repo_root, 'models'));
addpath(fullfile(repo_root, 'controllers'));
addpath(fullfile(repo_root, 'plotting'));

close all;
try
    bdclose('all');
catch
end

% Run MATLAB baseline
if nargin < 1
    P = [];
end
main_run_sim("baseline", P);

% Run Simulink baseline
run_simulink(P);

% Compare + overlay
compare_matlab_vs_simulink();
end
