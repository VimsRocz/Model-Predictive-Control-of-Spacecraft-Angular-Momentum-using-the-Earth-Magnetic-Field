function plot_results_simulink(matfile)
%PLOT_RESULTS_SIMULINK Plot Simulink outputs saved from run_simulink.

S = load(matfile);
P = S.P;

x_log = S.x_log;
m_log = S.m_log;
B_log = S.B_log;
tau_ext_log = S.tau_ext_log;
tau_mtq_log = S.tau_mtq_log;
tau_total_log = S.tau_total_log;

% Extract time and values
t = x_log.time(:);
x = x_log.signals.values;
m = m_log.signals.values;
B = B_log.signals.values;
tau_ext = tau_ext_log.signals.values;
tau_mtq = tau_mtq_log.signals.values;
tau_total = tau_total_log.signals.values;

fig1 = plot_time_series_report(P, t, x, m, B, tau_ext, tau_mtq, tau_total, "Simulink");
fig2 = plot_earth_field_scene(P, t, x, m, B, "Simulink");

out_dir = fileparts(matfile);
exportgraphics(fig1, fullfile(out_dir, "report_timeseries_Simulink.png"), "Resolution", 250);
exportgraphics(fig2, fullfile(out_dir, "scene_earth_Simulink.png"), "Resolution", 250);
end
