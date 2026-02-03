clear; clc;

P = params_default();

% Choose controller: "baseline" or "mpc"
ctrlMode = "mpc";

x = zeros(3, P.N);     % momentum error state
m = zeros(3, P.N);     % dipole command
tau_ext_log = zeros(3, P.N);
B_log = zeros(3, P.N);

% Initial momentum (example)
x(:,1) = [0.02; -0.01; 0.03]; % [N*m*s] ~ wheel momentum error

% Providers
B_provider   = @(kk) orbit_propagator_simple(kk, P);
tau_provider = @(kk) ext_torques(kk, P);

m_prev = zeros(3,1);

for k = 1:P.N-1
    Bk = B_provider(k);
    tk = tau_provider(k);

    B_log(:,k) = Bk;
    tau_ext_log(:,k) = tk;

    if ctrlMode == "baseline"
        mk = baseline_mtq(x(:,k), Bk, P);
    else
        mk = mpc_mtq_qp(x(:,k), k, P, B_provider, tau_provider, m_prev);
    end

    m(:,k) = mk;

    % Plant update: x_{k+1} = x_k + Ts*(tau_ext - [B]_x*m)
    x(:,k+1) = x(:,k) + P.Ts*( tk - skew(Bk)*mk );

    m_prev = mk;
end

save("sim_out.mat","P","x","m","tau_ext_log","B_log","ctrlMode");

plot_results("sim_out.mat");
