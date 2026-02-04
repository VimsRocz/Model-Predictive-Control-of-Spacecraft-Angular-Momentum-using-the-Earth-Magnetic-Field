function compare_matlab_vs_simulink()
%COMPARE_MATLAB_VS_SIMULINK Compare MATLAB and Simulink outputs step-by-step.

repo_root = fileparts(mfilename('fullpath'));
matlab_file = fullfile(repo_root, 'MATLAB', 'Output', 'momentum', 'sim_out_matlab.mat');
simulink_file = fullfile(repo_root, 'Simulink', 'Output', 'momentum', 'sim_out_simulink.mat');

S1 = load(matlab_file);
S2 = load(simulink_file);

P = S1.P;

fprintf('MATLAB file:   %s\n', matlab_file);
fprintf('Simulink file: %s\n\n', simulink_file);

% MATLAB timeline
Tm = (0:P.N-1)*P.Ts;
Xm = S1.x.'; % Nx3
Mm = S1.m.'; % Nx3

% Simulink logs
x_log = S2.x_log; m_log = S2.m_log;
Ts = x_log.time; Xs = x_log.signals.values; Ms = m_log.signals.values;
if size(Xs,2) ~= 3 && size(Xs,1) == 3
    Xs = Xs.';
end
if size(Ms,2) ~= 3 && size(Ms,1) == 3
    Ms = Ms.';
end

% Interpolate Simulink data onto MATLAB time grid
Xi = interp1(Ts, Xs, Tm, 'linear', 'extrap');
Mi = interp1(Ts, Ms, Tm, 'linear', 'extrap');

x_diff = vecnorm(Xm - Xi, 2, 2);
m_diff = vecnorm(Mm - Mi, 2, 2);

[dx_max, ix] = max(x_diff);
[dm_max, im] = max(m_diff);

fprintf('Max |x_diff| = %.3e at t=%.1f s\n', dx_max, Tm(ix));
fprintf('Max |m_diff| = %.3e at t=%.1f s\n', dm_max, Tm(im));

fprintf('\nInitial conditions:\n');
fprintf('  MATLAB   x(0) = [% .4e % .4e % .4e]\n', Xm(1,1), Xm(1,2), Xm(1,3));
fprintf('  Simulink x(0) = [% .4e % .4e % .4e]\n', Xi(1,1), Xi(1,2), Xi(1,3));
fprintf('  MATLAB   m(0) = [% .4e % .4e % .4e]\n', Mm(1,1), Mm(1,2), Mm(1,3));
fprintf('  Simulink m(0) = [% .4e % .4e % .4e]\n\n', Mi(1,1), Mi(1,2), Mi(1,3));

% Find first time diff exceeds tolerance
xtol = 1e-6; mtol = 1e-6;
ix1 = find(x_diff > xtol, 1, 'first');
im1 = find(m_diff > mtol, 1, 'first');

if isempty(ix1)
    fprintf('x_diff never exceeds %.1e\n', xtol);
else
    fprintf('x_diff first exceeds %.1e at t=%.1f s, |x|=%.3e\n', xtol, Tm(ix1), x_diff(ix1));
end
if isempty(im1)
    fprintf('m_diff never exceeds %.1e\n', mtol);
else
    fprintf('m_diff first exceeds %.1e at t=%.1f s, |m|=%.3e\n', mtol, Tm(im1), m_diff(im1));
end

% Quick overlay plots
figure; plot(Tm, vecnorm(Xm,2,2), 'b', Tm, vecnorm(Xi,2,2), 'r--'); grid on;
legend('MATLAB','Simulink');
xlabel('time [s]'); ylabel('||x|| [N·m·s]');
title('Momentum magnitude comparison');

figure; plot(Tm, Mm, 'LineWidth', 1.0); hold on; plot(Tm, Mi, '--', 'LineWidth', 1.0); grid on;
legend('m_x (MAT)','m_y (MAT)','m_z (MAT)','m_x (SIM)','m_y (SIM)','m_z (SIM)');
xlabel('time [s]'); ylabel('m [A·m^2]');
title('Dipole command comparison');
end
