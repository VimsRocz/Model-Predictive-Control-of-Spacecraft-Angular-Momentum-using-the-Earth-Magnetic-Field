function interactive_mtq_explorer()
%INTERACTIVE_MTQ_EXPLORER Simple UI to vary orbit/MTQ parameters and rerun simulations.
%
% Buttons:
%   - Run MATLAB (baseline)
%   - Run Simulink (baseline)
%   - Run Both + Compare
%   - Animate last MATLAB / Simulink run
%
% Note: The Simulink model is rebuilt on each run to embed parameters (keeps results reproducible).

repo_root = fileparts(mfilename('fullpath'));
addpath(repo_root);
addpath(fullfile(repo_root, 'params'));
addpath(fullfile(repo_root, 'models'));
addpath(fullfile(repo_root, 'controllers'));
addpath(fullfile(repo_root, 'plotting'));

P0 = params_default();

fig = uifigure('Name','MTQ Explorer (Baseline vs MPC + Visualization)', 'Position',[100 100 620 920]);
gl = uigridlayout(fig, [32 2]);
gl.ColumnWidth = {190,'1x'};
gl.RowHeight = repmat({32}, 1, 31);
gl.RowHeight{32} = '1x';
gl.Padding = [12 12 12 12];

% --- Parameter fields ---
uilabel(gl, 'Text','Plant model', 'HorizontalAlignment','left');
plantMode = uidropdown(gl, 'Items',{'momentum-only','full (attitude+wheels)'}, 'Value','momentum-only');

uilabel(gl, 'Text','Controller (MATLAB)', 'HorizontalAlignment','left');
ctrlMode = uidropdown(gl, 'Items',{'baseline','mpc'}, 'Value','baseline');

uilabel(gl, 'Text','B-field model', 'HorizontalAlignment','left');
bfieldModel = uidropdown(gl, 'Items',{'igrf','dipole'}, 'Value', char(P0.env.bfield_model));

uilabel(gl, 'Text','IGRF decimal year', 'HorizontalAlignment','left');
igrfYear = uieditfield(gl, 'numeric', 'Value', P0.env.igrf_decimal_year, 'Limits',[1900 2100]);

uilabel(gl, 'Text','Orbit altitude [km]', 'HorizontalAlignment','left');
alt_km = uieditfield(gl, 'numeric', 'Value', P0.orbit.alt_m/1e3, 'Limits',[100 4000], 'LowerLimitInclusive','on');

uilabel(gl, 'Text','Inclination [deg]', 'HorizontalAlignment','left');
inc_deg = uieditfield(gl, 'numeric', 'Value', P0.orbit.inc_deg, 'Limits',[0 180], 'LowerLimitInclusive','on');

uilabel(gl, 'Text','Initial anomaly θ0 [deg]', 'HorizontalAlignment','left');
theta0_deg = uieditfield(gl, 'numeric', 'Value', P0.orbit.theta0_deg, 'Limits',[-360 360]);

uilabel(gl, 'Text','B0 (equator surface) [µT]', 'HorizontalAlignment','left');
B0_uT = uieditfield(gl, 'numeric', 'Value', 1e6*P0.B0_T, 'Limits',[1 200]);

uilabel(gl, 'Text','m_max (per axis) [A·m^2]', 'HorizontalAlignment','left');
mmax = uieditfield(gl, 'numeric', 'Value', P0.m_max(1), 'Limits',[0.001 5]);

uilabel(gl, 'Text','KH gain (diag) [1/s]', 'HorizontalAlignment','left');
KH_gain = uieditfield(gl, 'numeric', 'Value', P0.KH(1,1), 'Limits',[1e-6 1e-1]);

uilabel(gl, 'Text','MPC horizon Nh [steps]', 'HorizontalAlignment','left');
mpcNh = uieditfield(gl, 'numeric', 'Value', P0.mpc.Nh, 'Limits',[5 300]);

uilabel(gl, 'Text','Enable bias torque', 'HorizontalAlignment','left');
chkBias = uicheckbox(gl, 'Value', logical(P0.env.enable_bias));

uilabel(gl, 'Text','Enable periodic torque', 'HorizontalAlignment','left');
chkPeriodic = uicheckbox(gl, 'Value', logical(P0.env.enable_periodic));

uilabel(gl, 'Text','Enable gravity gradient', 'HorizontalAlignment','left');
chkGG = uicheckbox(gl, 'Value', false);

uilabel(gl, 'Text','GG scale factor', 'HorizontalAlignment','left');
ggScale = uieditfield(gl, 'numeric', 'Value', P0.env.gg_scale, 'Limits',[0 1e6]);

uilabel(gl, 'Text','Enable aerodynamic drag', 'HorizontalAlignment','left');
chkDrag = uicheckbox(gl, 'Value', false);

uilabel(gl, 'Text','Drag scale factor', 'HorizontalAlignment','left');
dragScale = uieditfield(gl, 'numeric', 'Value', P0.env.drag_scale, 'Limits',[0 1e6]);

uilabel(gl, 'Text','Initial x0_x [N·m·s]', 'HorizontalAlignment','left');
x0x = uieditfield(gl, 'numeric', 'Value', P0.x0(1));

uilabel(gl, 'Text','Initial x0_y [N·m·s]', 'HorizontalAlignment','left');
x0y = uieditfield(gl, 'numeric', 'Value', P0.x0(2));

uilabel(gl, 'Text','Initial x0_z [N·m·s]', 'HorizontalAlignment','left');
x0z = uieditfield(gl, 'numeric', 'Value', P0.x0(3));

uilabel(gl, 'Text','ADCS Kp (full)', 'HorizontalAlignment','left');
adcsKp = uieditfield(gl, 'numeric', 'Value', P0.adcs.Kp, 'Limits',[0 10]);

uilabel(gl, 'Text','ADCS Kd (full)', 'HorizontalAlignment','left');
adcsKd = uieditfield(gl, 'numeric', 'Value', P0.adcs.Kd, 'Limits',[0 10]);

uilabel(gl, 'Text','RW torque max [N·m] (full)', 'HorizontalAlignment','left');
rwTauMax = uieditfield(gl, 'numeric', 'Value', P0.rw.tau_max_Nm(1), 'Limits',[1e-6 1]);

uilabel(gl, 'Text','Full int substeps', 'HorizontalAlignment','left');
intSub = uieditfield(gl, 'numeric', 'Value', P0.full.int_substeps, 'Limits',[1 100]);

uilabel(gl, 'Text','Orbit markers (3D) [count]', 'HorizontalAlignment','left');
nMarkers = uieditfield(gl, 'numeric', 'Value', 36, 'Limits',[6 200]);

% --- Buttons ---
btnRunMatlab = uibutton(gl, 'Text','Run MATLAB (selected)', 'ButtonPushedFcn', @onRunMatlab);
btnRunMatlab.Layout.Row = 26; btnRunMatlab.Layout.Column = 1;

btnRunSimulink = uibutton(gl, 'Text','Run Simulink', 'ButtonPushedFcn', @onRunSimulink);
btnRunSimulink.Layout.Row = 26; btnRunSimulink.Layout.Column = 2;

btnRunBoth = uibutton(gl, 'Text','MATLAB vs Simulink (baseline)', 'ButtonPushedFcn', @onRunBoth);
btnRunBoth.Layout.Row = 27; btnRunBoth.Layout.Column = 1;

btnRunMatlabBoth = uibutton(gl, 'Text','Baseline vs MPC (MATLAB)', 'ButtonPushedFcn', @onRunMatlabBothControllers);
btnRunMatlabBoth.Layout.Row = 27; btnRunMatlabBoth.Layout.Column = 2;

btnAnimMatlab = uibutton(gl, 'Text','Animate MATLAB (selected)', 'ButtonPushedFcn', @onAnimateMatlab);
btnAnimMatlab.Layout.Row = 28; btnAnimMatlab.Layout.Column = 1;

btnAnimSim = uibutton(gl, 'Text','Animate Simulink', 'ButtonPushedFcn', @onAnimateSimulink);
btnAnimSim.Layout.Row = 28; btnAnimSim.Layout.Column = 2;

btnOpenOut = uibutton(gl, 'Text','Open outputs/', 'ButtonPushedFcn', @onOpenOutputs);
btnOpenOut.Layout.Row = 29; btnOpenOut.Layout.Column = 1;

btnBench = uibutton(gl, 'Text','Run Benchmarks', 'ButtonPushedFcn', @onRunBenchmarks);
btnBench.Layout.Row = 29; btnBench.Layout.Column = 2;

btnAnimFull = uibutton(gl, 'Text','Animate FULL (last)', 'ButtonPushedFcn', @onAnimateFull);
btnAnimFull.Layout.Row = 30; btnAnimFull.Layout.Column = [1 2];

% --- Status log ---
status = uitextarea(gl, 'Editable','off');
status.Layout.Row = 32;
status.Layout.Column = [1 2];
status.Value = { ...
    'Adjust parameters, then run MATLAB/Simulink.' ...
    'Tip: "Baseline vs MPC (MATLAB)" runs both controllers with the same environment.' ...
    'Tip: choose plant model "full" to simulate LVLH attitude + wheels + MTQ desaturation.' ...
    '' ...
    };

% ---------- helpers ----------
function logLine(msg)
    status.Value = [status.Value; {char(string(msg))}]; %#ok<AGROW>
    drawnow limitrate;
end

function P = readP()
    P = params_default(); % start from known-good defaults
    P.env.bfield_model = string(bfieldModel.Value);
    P.env.igrf_decimal_year = igrfYear.Value;
    P.orbit.alt_m = alt_km.Value * 1e3;
    P.orbit.inc_deg = inc_deg.Value;
    P.orbit.theta0_deg = theta0_deg.Value;
    P.B0_T = B0_uT.Value * 1e-6;
    P.m_max = mmax.Value * ones(3,1);
    P.KH = KH_gain.Value * eye(3);
    P.mpc.Nh = max(5, round(mpcNh.Value));
    P.env.enable_bias = logical(chkBias.Value);
    P.env.enable_periodic = logical(chkPeriodic.Value);
    P.env.enable_gravity_gradient = logical(chkGG.Value);
    P.env.gg_scale = ggScale.Value;
    P.env.enable_drag = logical(chkDrag.Value);
    P.env.drag_scale = dragScale.Value;
    P.x0 = [x0x.Value; x0y.Value; x0z.Value];
    P.adcs.Kp = adcsKp.Value;
    P.adcs.Kd = adcsKd.Value;
    P.rw.tau_max_Nm = rwTauMax.Value * ones(3,1);
    P.full.int_substeps = max(1, round(intSub.Value));
    P.N = round(P.Tend / P.Ts);
    P.viz.nOrbitMarkers = nMarkers.Value;
end

function out = matlab_out_file(mode)
    if nargin < 1 || strlength(mode) == 0
        mode = string(ctrlMode.Value);
    end
    safe = string(regexprep(char(mode), '[^A-Za-z0-9]+', '_'));
    out = fullfile(repo_root, 'outputs', 'matlab', "sim_out_matlab_" + safe + ".mat");
    if ~exist(out, 'file')
        out = fullfile(repo_root, 'outputs', 'matlab', 'sim_out_matlab.mat');
    end
end

function out = simulink_out_file()
    out = fullfile(repo_root, 'outputs', 'simulink', 'sim_out_simulink.mat');
end

function out = simulink_full_out_file()
    out = fullfile(repo_root, 'outputs', 'simulink_full', 'sim_out_full_simulink.mat');
end

% ---------- callbacks ----------
function onRunMatlab(~, ~)
    P = readP();
    mode = string(ctrlMode.Value);
    logLine("Running MATLAB (" + mode + ", " + string(plantMode.Value) + ")…");
    try
        if startsWith(string(plantMode.Value), "full")
            main_run_full_sim(mode, P);
            logLine("MATLAB full complete: " + full_out_file(mode));
        else
            main_run_sim(mode, P);
            logLine("MATLAB complete: " + matlab_out_file(mode));
        end
    catch ME
        logLine("MATLAB error: " + ME.message);
        uialert(fig, ME.message, 'MATLAB Run Failed');
    end
end

function onRunSimulink(~, ~)
    P = readP();
    logLine("Running Simulink (" + string(plantMode.Value) + ")…");
    try
        if startsWith(string(plantMode.Value), "full")
            run_simulink_full(P);
            logLine("Simulink full complete: " + simulink_full_out_file());
        else
            run_simulink(P);
            logLine("Simulink complete: " + simulink_out_file());
        end
    catch ME
        logLine("Simulink error: " + ME.message);
        uialert(fig, ME.message, 'Simulink Run Failed');
    end
end

function onRunBoth(~, ~)
    P = readP();
    logLine("Running MATLAB + Simulink (" + string(plantMode.Value) + ")…");
    try
        if startsWith(string(plantMode.Value), "full")
            main_run_full_sim("baseline", P);
            run_simulink_full(P);
            logLine("Full MATLAB + Simulink complete. See outputs/matlab_full and outputs/simulink_full.");
        else
            run_all_and_compare(P);
            logLine("Compare complete. See outputs/ and console for diffs.");
        end
    catch ME
        logLine("Run/compare error: " + ME.message);
        uialert(fig, ME.message, 'Run Both Failed');
    end
end

function onRunMatlabBothControllers(~, ~)
    P = readP();
    logLine("Running MATLAB baseline + MPC (" + string(plantMode.Value) + ")…");
    try
        if startsWith(string(plantMode.Value), "full")
            main_run_full_sim("baseline", P);
            main_run_full_sim("mpc", P);
            logLine("Baseline vs MPC full complete. See outputs/matlab_full/ for files.");
        else
            main_run_sim("baseline", P);
            main_run_sim("mpc", P);
            compare_baseline_vs_mpc(P);
            logLine("Baseline vs MPC complete. See outputs/matlab/ for files.");
        end
    catch ME
        logLine("Run/compare error: " + ME.message);
        uialert(fig, ME.message, 'Baseline vs MPC Failed');
    end
end

function onAnimateMatlab(~, ~)
    if startsWith(string(plantMode.Value), "full")
        onAnimateFull();
        return;
    end

    f = matlab_out_file(string(ctrlMode.Value));
    if ~exist(f, 'file')
        uialert(fig, "MATLAB output not found. Run MATLAB first.", 'Missing Output');
        return;
    end
    S = load(f);
    P = S.P;
    t = (0:size(S.x,2)-1) * P.Ts;
    lbl = "MATLAB/baseline";
    if isfield(S,'ctrlMode')
        lbl = "MATLAB/" + string(S.ctrlMode);
    end
    animate_earth_field_scene(P, t, S.x.', S.m.', S.B_log.', lbl);
end

function out = full_out_file(mode)
    if nargin < 1 || strlength(mode) == 0
        mode = string(ctrlMode.Value);
    end
    safe = string(regexprep(char(mode), '[^A-Za-z0-9]+', '_'));
    out = fullfile(repo_root, 'outputs', 'matlab_full', "sim_out_full_" + safe + ".mat");
    if ~exist(out, 'file')
        out = fullfile(repo_root, 'outputs', 'matlab_full', 'sim_out_full.mat');
    end
end

function onAnimateFull(~, ~)
    f = full_out_file(string(ctrlMode.Value));
    if ~exist(f, 'file')
        uialert(fig, "Full MATLAB output not found. Run MATLAB with plant model 'full' first.", 'Missing Output');
        return;
    end
    S = load(f);
    P = S.P;
    animate_mtq_full_scene(P, S, "MATLAB/full/" + string(S.desatMode));
end

function onAnimateSimulink(~, ~)
    if startsWith(string(plantMode.Value), "full")
        f = simulink_full_out_file();
        if ~exist(f, 'file')
            uialert(fig, "Simulink full output not found. Run Simulink (full) first.", 'Missing Output');
            return;
        end
        S = load(f);
        P = S.P;
        animate_mtq_full_scene(P, S, "Simulink/full");
        return;
    end

    f = simulink_out_file();
    if ~exist(f, 'file')
        uialert(fig, "Simulink output not found. Run Simulink first.", 'Missing Output');
        return;
    end
    S = load(f);
    P = S.P;
    x = S.x_log.signals.values;
    m = S.m_log.signals.values;
    B = S.B_log.signals.values;
    t = S.x_log.time;
    animate_earth_field_scene(P, t, x, m, B, "Simulink");
end

function onOpenOutputs(~, ~)
    out_dir = fullfile(repo_root, 'outputs');
    if exist(out_dir, 'dir')
        try
            open(out_dir);
        catch
            logLine("Outputs dir: " + out_dir);
        end
    end
end

function onRunBenchmarks(~, ~)
    logLine("Running thesis-style scenario benchmarks (MATLAB)…");
    try
        run_thesis_benchmarks();
        logLine("Benchmarks complete. See outputs/analysis/.");
    catch ME
        logLine("Benchmarks error: " + ME.message);
        uialert(fig, ME.message, 'Benchmarks Failed');
    end
end
end
