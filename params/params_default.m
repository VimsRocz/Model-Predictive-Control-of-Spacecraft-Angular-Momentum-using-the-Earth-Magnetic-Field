function P = params_default()
%PARAMS_DEFAULT Central parameter struct

% --- Simulation ---
P.Ts   = 1.0;           % [s] controller sample time
P.Tend = 2*3600;        % [s] total simulation time
P.N    = round(P.Tend/P.Ts);
P.x0   = [0.02; -0.01; 0.03]; % [N*m*s] initial momentum error

% --- Earth constants (used by simple orbit + dipole B-field visualizations) ---
P.earth.Re_m = 6371e3;              % [m] mean Earth radius
P.earth.mu_m3s2 = 3.986004418e14;   % [m^3/s^2] Earth gravitational parameter
P.earth.omega_radps = 7.2921159e-5; % [rad/s] Earth rotation rate
P.earth.dipole_tilt_deg = 11.0;     % [deg] simple dipole tilt (not IGRF)

% --- MTQ limits ---
P.m_max = [0.2; 0.2; 0.2];   % [A*m^2] per-axis dipole saturation (example)

% --- Spacecraft properties (used by external torque models) ---
P.sc.I_kgm2 = diag([0.02 0.018 0.015]); % [kg*m^2] inertia matrix (example CubeSat)
P.sc.Cd = 2.2;                           % [-] drag coefficient
P.sc.area_m2 = 0.01;                     % [m^2] reference area
P.sc.cp_offset_m = [0.02; 0; 0];         % [m] center-of-pressure offset (body/inertial, simplified)

% --- Environment / external torques ---
P.tau_ext_const = [2e-6; -1e-6; 1e-6]; % [N*m] constant bias torque example
P.env.enable_bias = true;
P.env.enable_periodic = true;
P.env.enable_gravity_gradient = false;
P.env.enable_drag = false;

% --- Earth magnetic field model ---
% "igrf" uses Aerospace Toolbox igrfmagm (IGRF-14), "dipole" uses a simple tilted dipole.
P.env.bfield_model = "igrf";
P.env.igrf_generation = 14;
P.env.igrf_decimal_year = 2026.0; % constant over the short simulation (hours)

% Drag density model: rho = rho_ref * exp(-(alt-alt_ref)/H)
P.env.drag_rho_ref_kgm3 = 1e-12;    % [kg/m^3] reference density (order-of-magnitude, LEO)
P.env.drag_alt_ref_m = 400e3;       % [m]
P.env.drag_scale_height_m = 55e3;   % [m]
P.env.drag_scale = 1.0;             % [-] scale factor for "what-if" sweeps

% Gravity-gradient scale factor (what-if sweeps)
P.env.gg_scale = 1.0;

% --- Orbit / B-field ---
P.orbit.alt_m  = 550e3;     % [m]
P.orbit.inc_deg = 97.5;     % [deg]
P.orbit.lat0_deg = 0;       % [deg] start latitude approx
P.B0_T = 30e-6;             % [T] rough magnitude (use real model if possible)
P.orbit.theta0_deg = 0;     % [deg] initial true anomaly (used in simple orbit)

% --- Visualization defaults ---
P.viz.nOrbitMarkers = 36;
P.viz.sat_size_Re = 0.02; % spacecraft cube size as fraction of Earth radius (for visibility)

% --- Attitude control + reaction wheel limits (full simulation) ---
P.adcs.Kp = 0.05;              % PD attitude gain (unitless)
P.adcs.Kd = 0.2;               % PD rate gain [N*m*s]
P.rw.tau_max_Nm = 5e-4 * ones(3,1); % [N*m] per-axis wheel torque limit (example)

% --- Simulink defaults ---
P.simulink.controller = "baseline"; % "baseline" or "mpc" for interpreted Simulink model
P.full.int_substeps = 10;           % internal integration substeps for full simulation

% --- Baseline controller gain ---
P.KH = 1e-3 * eye(3);       % momentum reduction gain

% --- MPC ---
P.mpc.Nh = 300;             % horizon steps (long enough to see B-field geometry change)
P.mpc.Q  = 1e4 * eye(3);    % emphasize momentum reduction
P.mpc.R  = 1e-2 * eye(3);   % small control-usage penalty (MTQs are cheap vs momentum growth)
P.mpc.S  = 1e-2 * eye(3);   % delta-m penalty (smoothness)
P.mpc.useDelta = true;
P.mpc.update_every = 20;     % [steps] solve QP every N steps, hold command between solves
P.mpc.warm_start = true;     % reuse previous solution as quadprog initial guess
P.mpc.enforce_tau_perpB = false; % optional: reduce m-component along B without changing tau (see mpc_mtq_qp)

% Solver options (quadprog)
P.mpc.qp_opts = optimoptions('quadprog', ...
    'Display','off', ...
    'Algorithm','interior-point-convex');

end
