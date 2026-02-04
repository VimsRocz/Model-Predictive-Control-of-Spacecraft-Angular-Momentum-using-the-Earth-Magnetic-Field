function tau_ext = ext_torques(k, P)
%EXT_TORQUES External torque model (bias + optional periodic + GG + drag).
%
% This is a simplified environment model to support controller benchmarking.
% All vectors are expressed in the same inertial/ECI frame (body=inertial assumption).

t_s = (k-1)*P.Ts;

tau_ext = zeros(3,1);

% --- Bias torque (e.g., residual dipole, offsets) ---
if ~isfield(P,'env') || ~isfield(P.env,'enable_bias') || P.env.enable_bias
    tau_ext = tau_ext + P.tau_ext_const(:);
end

% --- Small periodic content (proxy for varying aero/illumination) ---
if isfield(P,'env') && isfield(P.env,'enable_periodic') && P.env.enable_periodic
    tau_ext = tau_ext + 1e-6*[sin(2*pi*t_s/600); cos(2*pi*t_s/800); sin(2*pi*t_s/1000)];
end

% Orbit state (for GG + drag)
[r_eci, v_eci] = orbit_state_rv_simple(k, P);
rn = norm(r_eci);
Re = P.earth.Re_m;
alt = max(0, rn - Re);

% --- Gravity gradient torque: tau = 3*mu/r^3 * (rhat x (I*rhat)) ---
if isfield(P,'env') && isfield(P.env,'enable_gravity_gradient') && P.env.enable_gravity_gradient
    mu = P.earth.mu_m3s2;
    I = P.sc.I_kgm2;
    rhat = r_eci / max(rn, 1.0);
    tau_gg = 3 * mu / max(rn^3, 1.0) * cross(rhat, I*rhat);
    gg_scale = 1.0;
    if isfield(P.env,'gg_scale'); gg_scale = P.env.gg_scale; end
    tau_ext = tau_ext + gg_scale * tau_gg;
end

% --- Aerodynamic drag torque: tau = r_cp x F_drag ---
if isfield(P,'env') && isfield(P.env,'enable_drag') && P.env.enable_drag
    rho_ref = P.env.drag_rho_ref_kgm3;
    alt_ref = P.env.drag_alt_ref_m;
    H = P.env.drag_scale_height_m;
    drag_scale = 1.0;
    if isfield(P.env,'drag_scale'); drag_scale = P.env.drag_scale; end

    rho = rho_ref * exp(-(alt - alt_ref)/max(H, 1.0));
    rho = max(rho, 0);

    v = v_eci;
    vmag = norm(v);
    if vmag > 1e-9
        vhat = v / vmag;
        Fdrag = -0.5 * rho * P.sc.Cd * P.sc.area_m2 * vmag^2 * vhat;
        tau_drag = cross(P.sc.cp_offset_m(:), Fdrag);
        tau_ext = tau_ext + drag_scale * tau_drag;
    end
end
end
