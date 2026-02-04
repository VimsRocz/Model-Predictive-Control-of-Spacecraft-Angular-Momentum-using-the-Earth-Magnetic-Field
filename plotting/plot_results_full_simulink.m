function plot_results_full_simulink(matfile)
%PLOT_RESULTS_FULL_SIMULINK Plot full-plant Simulink outputs.

S = load(matfile);
P = S.P;

label = "Simulink/full";

fig1 = plot_time_series_report_full(P, S, label);
fig2 = plot_earth_field_scene_full(P, S, label);

out_dir = fileparts(matfile);
exportgraphics(fig1, fullfile(out_dir, "report_full_timeseries_Simulink.png"), "Resolution", 250);
exportgraphics(fig2, fullfile(out_dir, "scene_full_earth_Simulink.png"), "Resolution", 250);
end

