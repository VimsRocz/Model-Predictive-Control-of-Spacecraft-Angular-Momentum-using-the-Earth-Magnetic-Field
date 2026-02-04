function plot_results_full(matfile)
%PLOT_RESULTS_FULL Load a full-simulation MAT file and generate plots.

S = load(matfile);

if ~isfield(S, 'P')
    error("MAT file missing P.");
end
P = S.P;

label = "MATLAB/full";
if isfield(S,'desatMode')
    label = "MATLAB/full/" + string(S.desatMode);
end

fig1 = plot_time_series_report_full(P, S, label);
fig2 = plot_earth_field_scene_full(P, S, label);

out_dir = fileparts(matfile);
safe_label = string(regexprep(char(label), '[^A-Za-z0-9]+', '_'));
exportgraphics(fig1, fullfile(out_dir, "report_full_timeseries_" + safe_label + ".png"), "Resolution", 250);
exportgraphics(fig2, fullfile(out_dir, "scene_full_earth_" + safe_label + ".png"), "Resolution", 250);
end

