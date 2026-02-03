function P = params_default()
%PARAMS_DEFAULT Central parameter struct

% --- Simulation ---
P.Ts   = 1.0;           % [s] controller sample time
P.Tend = 2*3600;        % [s] total simulation time
P.N    = round(P.Tend/P.Ts);

% --- MTQ limits ---
P.m_max = [0.2; 0.2; 0.2];   % [A*m^2] per-axis dipole saturation (example)

% --- Environment (simple) ---
P.tau_ext_const = [2e-6; -1e-6; 1e-6]; % [N*m] constant bias torque example

% --- Orbit / B-field ---
P.orbit.alt_m  = 550e3;     % [m]
P.orbit.inc_deg = 97.5;     % [deg]
P.orbit.lat0_deg = 0;       % [deg] start latitude approx
P.B0_T = 30e-6;             % [T] rough magnitude (use real model if possible)

% --- Baseline controller gain ---
P.KH = 1e-3 * eye(3);       % momentum reduction gain

% --- MPC ---
P.mpc.Nh = 60;              % horizon steps
P.mpc.Q  = 1e2 * eye(3);
P.mpc.R  = 1e0 * eye(3);
P.mpc.S  = 1e0 * eye(3);    % delta-m penalty
P.mpc.useDelta = true;

% Solver options (quadprog)
P.mpc.qp_opts = optimoptions('quadprog', ...
    'Display','off', ...
    'Algorithm','interior-point-convex');

end
