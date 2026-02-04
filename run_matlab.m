function run_matlab(P)
%RUN_MATLAB Single entry-point for MATLAB simulation (momentum-only or full).
%
% Uses:
%   P.plant.model = "momentum" | "full"
%   P.controller  = "baseline" | "mpc"
%
% Default: full model + baseline controller (see params_default.m)

clc;

repo_root = fileparts(mfilename('fullpath'));
matlab_proc = fullfile(repo_root, 'MATLAB', 'Processing');
addpath(matlab_proc);
addpath(fullfile(repo_root, 'MATLAB', 'Input', 'params'));
addpath(fullfile(matlab_proc, 'models'));
addpath(fullfile(matlab_proc, 'controllers'));
addpath(fullfile(matlab_proc, 'plotting'));

if nargin < 1 || isempty(P)
    P = params_default();
end

if ~isfield(P, 'plant') || ~isfield(P.plant, 'model')
    P.plant.model = "full";
end
if ~isfield(P, 'controller') || strlength(P.controller) == 0
    P.controller = "baseline";
end

P.simulink.controller = string(P.controller); % keep MATLAB/Simulink aligned

plant = lower(string(P.plant.model));
ctrl = lower(string(P.controller));

switch plant
    case "full"
        main_run_full_sim(ctrl, P);
        if isfield(P,'viz') && isfield(P.viz,'auto_animate_full') && P.viz.auto_animate_full
            try
                f = fullfile(repo_root, 'MATLAB', 'Output', 'full', 'sim_out_full.mat');
                if exist(f,'file')
                    S = load(f);
                    animate_mtq_full_scene(P, S, "MATLAB/full/" + string(ctrl));
                end
            catch
            end
        end
    case "momentum"
        main_run_sim(ctrl, P);
        if isfield(P,'viz') && isfield(P.viz,'auto_animate_momentum') && P.viz.auto_animate_momentum
            try
                f = fullfile(repo_root, 'MATLAB', 'Output', 'momentum', 'sim_out_matlab.mat');
                if exist(f,'file')
                    S = load(f);
                    t = (0:size(S.x,2)-1) * S.P.Ts;
                    animate_earth_field_scene(S.P, t, S.x.', S.m.', S.B_log.', "MATLAB/momentum/" + string(ctrl));
                end
            catch
            end
        end
    otherwise
        error("Unknown P.plant.model '%s'. Use 'full' or 'momentum'.", plant);
end

% If a Simulink run exists, print a quick consistency check summary.
try
    compare_last_runs(P);
catch
end
end

function compare_last_runs(P)
% Compare last MATLAB vs Simulink outputs (if available).
repo_root = fileparts(mfilename('fullpath'));

if strcmpi(string(P.plant.model), "full")
    fmat = fullfile(repo_root, 'MATLAB', 'Output', 'full', 'sim_out_full.mat');
    fsim = fullfile(repo_root, 'Simulink', 'Output', 'full', 'sim_out_full_simulink.mat');
    if ~exist(fmat,'file') || ~exist(fsim,'file')
        return;
    end
    Sm = load(fmat);
    Ss = load(fsim);
    t = Sm.t(:);
    hm = Sm.h_w;
    hs = Ss.h_w;
    ts = Ss.t(:);
    hs_i = interp1(ts, hs, t, 'linear', 'extrap');
    diff = vecnorm(hm - hs_i, 2, 2);
    fprintf('MATLAB vs Simulink (full) max |h_w diff| = %.3e\n', max(diff));
else
    fmat = fullfile(repo_root, 'MATLAB', 'Output', 'momentum', 'sim_out_matlab.mat');
    fsim = fullfile(repo_root, 'Simulink', 'Output', 'momentum', 'sim_out_simulink.mat');
    if ~exist(fmat,'file') || ~exist(fsim,'file')
        return;
    end
    Sm = load(fmat);
    Ss = load(fsim);
    t = (0:Sm.P.N-1).' * Sm.P.Ts;
    xm = Sm.x.';
    xs = Ss.x_log.signals.values;
    ts = Ss.x_log.time;
    xs_i = interp1(ts, xs, t, 'linear', 'extrap');
    diff = vecnorm(xm - xs_i, 2, 2);
    fprintf('MATLAB vs Simulink (momentum) max |x diff| = %.3e\n', max(diff));
end
end
