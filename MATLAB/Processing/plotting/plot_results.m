function plot_results(matfile)
S = load(matfile);
P = S.P;

% Expect MATLAB simulation arrays as 3xN
x = S.x;
m = S.m;
B = S.B_log;
tau_ext = S.tau_ext_log;

if isfield(S, 'tau_mtq_log')
    tau_mtq = S.tau_mtq_log;
else
    tau_mtq = cross(m.', B.').';
end

if isfield(S, 'tau_total_log')
    tau_total = S.tau_total_log;
else
    tau_total = tau_ext + tau_mtq;
end

t = (0:size(x,2)-1) * P.Ts;

if isfield(S, 'ctrlMode')
    label = "MATLAB/" + string(S.ctrlMode);
else
    label = "MATLAB";
end

fig1 = plot_time_series_report(P, t, x.', m.', B.', tau_ext.', tau_mtq.', tau_total.', label);
fig2 = plot_earth_field_scene(P, t, x.', m.', B.', label);

out_dir = fileparts(matfile);
safe_label = string(regexprep(char(label), '[^A-Za-z0-9]+', '_'));
exportgraphics(fig1, fullfile(out_dir, "report_timeseries_" + safe_label + ".png"), "Resolution", 250);
exportgraphics(fig2, fullfile(out_dir, "scene_earth_" + safe_label + ".png"), "Resolution", 250);
end
