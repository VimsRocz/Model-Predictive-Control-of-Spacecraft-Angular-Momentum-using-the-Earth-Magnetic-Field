function tau_ext_body = ext_torques_body(t_s, r_eci_m, v_eci_ms, q_ib, P)
%EXT_TORQUES_BODY External torque model expressed in the BODY frame.
%
% This is a more physically consistent variant of ext_torques.m for the full
% attitude + wheel simulation.
%
% Components:
% - Constant bias torque (body-fixed) + optional periodic content (body-fixed)
% - Optional gravity-gradient torque (uses r expressed in body frame)
% - Optional aerodynamic drag torque (uses F_drag in body frame and cp offset)
%
% Inputs
%   t_s      : time [s]
%   r_eci_m  : position in ECI [m]
%   v_eci_ms : velocity in ECI [m/s]
%   q_ib     : quaternion inertial->body (1x4, scalar-first)
%   P        : params
%
% Output
%   tau_ext_body : [3x1] external torque in body [N*m]

tau_ext_body = zeros(3,1);

% DCM inertial->body
C_bi = quat2dcm(q_ib);

% --- Bias torque (body-fixed constant) ---
if ~isfield(P,'env') || ~isfield(P.env,'enable_bias') || P.env.enable_bias
    tau_ext_body = tau_ext_body + P.tau_ext_const(:);
end

% --- Small periodic content (proxy for varying aero/illumination) ---
if isfield(P,'env') && isfield(P.env,'enable_periodic') && P.env.enable_periodic
    tau_ext_body = tau_ext_body + 1e-6*[sin(2*pi*t_s/600); cos(2*pi*t_s/800); sin(2*pi*t_s/1000)];
end

% Orbit geometry
r = r_eci_m(:);
v = v_eci_ms(:);
rn = norm(r);
Re = P.earth.Re_m;
alt = max(0, rn - Re);

% --- Gravity gradient torque (body frame): 3*mu/r^3 * (rhat_b x (I*rhat_b)) ---
if isfield(P,'env') && isfield(P.env,'enable_gravity_gradient') && P.env.enable_gravity_gradient
    mu = P.earth.mu_m3s2;
    I = P.sc.I_kgm2;
    rhat_b = (C_bi * (r / max(rn, 1.0)));
    tau_gg = 3 * mu / max(rn^3, 1.0) * cross(rhat_b, I*rhat_b);
    gg_scale = 1.0;
    if isfield(P.env,'gg_scale'); gg_scale = P.env.gg_scale; end
    tau_ext_body = tau_ext_body + gg_scale * tau_gg;
end

% --- Aerodynamic drag torque (body frame): tau = r_cp x F_drag_body ---
if isfield(P,'env') && isfield(P.env,'enable_drag') && P.env.enable_drag
    rho_ref = P.env.drag_rho_ref_kgm3;
    alt_ref = P.env.drag_alt_ref_m;
    H = P.env.drag_scale_height_m;
    drag_scale = 1.0;
    if isfield(P.env,'drag_scale'); drag_scale = P.env.drag_scale; end

    rho = rho_ref * exp(-(alt - alt_ref)/max(H, 1.0));
    rho = max(rho, 0);

    vmag = norm(v);
    if vmag > 1e-9
        vhat = v / vmag;
        Fdrag_eci = -0.5 * rho * P.sc.Cd * P.sc.area_m2 * vmag^2 * vhat;
        Fdrag_body = C_bi * Fdrag_eci;
        tau_drag = cross(P.sc.cp_offset_m(:), Fdrag_body);
        tau_ext_body = tau_ext_body + drag_scale * tau_drag;
    end
end
end

