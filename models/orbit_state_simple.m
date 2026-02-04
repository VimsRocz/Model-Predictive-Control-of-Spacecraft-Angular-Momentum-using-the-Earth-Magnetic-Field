function [r_eci_m, theta_rad] = orbit_state_simple(k, P)
%ORBIT_STATE_SIMPLE Very simple circular orbit state for visualization and B-field.
%
% Returns inertial position r_eci_m [m] and true anomaly theta_rad [rad]
% at discrete step k.

k = max(1, round(k));
t = (k-1) * P.Ts;

Re = P.earth.Re_m;
mu = P.earth.mu_m3s2;

r_orb = Re + P.orbit.alt_m;
inc = deg2rad(P.orbit.inc_deg);

if isfield(P.orbit, 'theta0_deg')
    theta0 = deg2rad(P.orbit.theta0_deg);
else
    theta0 = 0;
end

% Mean motion for circular orbit
w = sqrt(mu / r_orb^3);
theta_rad = theta0 + w * t;

% Orbit in x-y plane, then incline about x-axis.
r_orb_plane = r_orb * [cos(theta_rad); sin(theta_rad); 0];
R1 = [1 0 0; 0 cos(inc) -sin(inc); 0 sin(inc) cos(inc)];
r_eci_m = R1 * r_orb_plane;
end
