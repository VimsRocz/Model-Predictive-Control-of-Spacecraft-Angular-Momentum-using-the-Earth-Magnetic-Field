function plot_results(matfile)
P = [];
S = load(matfile);

P = S.P; x = S.x; m = S.m;
t = (0:P.N-1)*P.Ts;

figure; plot(t, vecnorm(x,2,1)); grid on;
xlabel('time [s]'); ylabel('||x|| [N·m·s]');
title("Momentum magnitude ("+S.ctrlMode+")");

figure; plot(t, m'); grid on;
xlabel('time [s]'); ylabel('m [A·m^2]');
legend('m_x','m_y','m_z');
title("Dipole commands ("+S.ctrlMode+")");
end
