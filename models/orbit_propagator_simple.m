function B_body = orbit_propagator_simple(k, P)
%ORBIT_PROPAGATOR_SIMPLE Earth magnetic field along a simple circular orbit.
%
% Momentum-only model assumes body frame ~ inertial frame, so this returns
% B in the same frame. The underlying model is selected via:
%   P.env.bfield_model = "dipole" | "igrf"

r_eci_m = orbit_state_simple(k, P);
t_s = (k-1) * P.Ts;
B_eci_T = earth_mag_field_eci(r_eci_m, P, t_s);

B_body = B_eci_T;
end
