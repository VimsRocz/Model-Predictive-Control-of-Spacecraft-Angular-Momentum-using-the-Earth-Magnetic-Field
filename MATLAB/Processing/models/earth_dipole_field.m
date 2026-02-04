function B_eci_T = earth_dipole_field(r_eci_m, P, t_s)
%EARTH_DIPOLE_FIELD Idealized Earth dipole magnetic field with optional tilt + Earth rotation.
%
% B_eci_T is in Tesla.
% P.B0_T is treated as equatorial surface magnitude (at r = Re, z=0).
%
% Optional:
%   P.earth.omega_radps         Earth rotation rate [rad/s] (default 7.2921159e-5)
%   P.earth.dipole_tilt_deg     dipole tilt relative to +Z (default 0)
%
% Note: This is not IGRF; it's a simple dipole for algorithm + visualization testing.

if nargin < 3
    t_s = 0;
end

Re = P.earth.Re_m;
B0 = P.B0_T;

r = r_eci_m(:);
rn = norm(r);

if rn < 1.0
    B_eci_T = zeros(3,1);
    return;
end

omega = 7.2921159e-5; % [rad/s]
if isfield(P, 'earth') && isfield(P.earth, 'omega_radps')
    omega = P.earth.omega_radps;
end
tilt = 0; % [rad]
if isfield(P, 'earth') && isfield(P.earth, 'dipole_tilt_deg')
    tilt = deg2rad(P.earth.dipole_tilt_deg);
end

% ECI -> ECEF rotation (about +Z)
ct = cos(omega * t_s);
st = sin(omega * t_s);
R3_eci2ecef = [ct st 0; -st ct 0; 0 0 1];

r_ecef = R3_eci2ecef * r;
r_hat = r_ecef / rn;

% Dipole axis in ECEF (tilt in X-Z plane at t=0 for simplicity)
m_hat = [sin(tilt); 0; cos(tilt)];
mr = dot(m_hat, r_hat);

% Dipole field: B = B0 * (Re/r)^3 * (3 r_hat (m_hat·r_hat) - m_hat)
B_ecef_T = B0 * (Re/rn)^3 * (3 * r_hat * mr - m_hat);

% ECEF -> ECI
B_eci_T = R3_eci2ecef.' * B_ecef_T;
end
