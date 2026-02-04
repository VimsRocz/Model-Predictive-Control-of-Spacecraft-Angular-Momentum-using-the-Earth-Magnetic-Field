function B_eci_T = earth_igrf_field(r_eci_m, P, t_s)
%EARTH_IGRF_FIELD Earth magnetic field using IGRF-14 (Aerospace Toolbox).
%
% Inputs
%   r_eci_m : [3x1] position in ECI [m]
%   P       : params struct (expects P.earth.omega_radps, P.earth.Re_m)
%   t_s     : time since epoch [s]
%
% Output
%   B_eci_T : [3x1] magnetic field in ECI [Tesla]
%
% Implementation notes
% - Converts ECI -> ECEF via simple Earth rotation about +Z:
%     r_ecef = R3(omega*t) * r_eci
% - Uses ecef2lla for geodetic lat/lon/alt
% - Calls igrfmagm(height_km, lat_deg, lon_deg, decimalYear, generation)
% - Converts NED -> ECEF (ned2ecefv) then ECEF -> ECI (R3')
%
% This is intended for simulation/visualization and does not include full
% Earth orientation (precession/nutation/polar motion).

if nargin < 3
    t_s = 0;
end

if exist('igrfmagm', 'file') ~= 2
    error("igrfmagm not found. Install Aerospace Toolbox or set P.env.bfield_model='dipole'.");
end

r_eci_m = r_eci_m(:);
if numel(r_eci_m) ~= 3
    error("r_eci_m must be 3x1.");
end

omega = P.earth.omega_radps;
ct = cos(omega * t_s);
st = sin(omega * t_s);
R3 = [ct st 0; -st ct 0; 0 0 1]; % ECI -> ECEF (simple rotation)

r_ecef_m = R3 * r_eci_m;

% ECEF -> geodetic LLA (deg, deg, m)
lla = ecef2lla(r_ecef_m.');
lat_deg = lla(1);
lon_deg = lla(2);
alt_m = lla(3);

% igrfmagm wants height in km
height_km = alt_m / 1e3;

% Decimal year: treat model as constant over short runs unless user provides epoch.
if isfield(P, 'env') && isfield(P.env, 'igrf_decimal_year') && ~isempty(P.env.igrf_decimal_year)
    dy = double(P.env.igrf_decimal_year);
elseif isfield(P, 'env') && isfield(P.env, 'epoch_utc') && ~isempty(P.env.epoch_utc)
    try
        dy = double(decyear(P.env.epoch_utc + seconds(t_s)));
    catch
        dy = double(decyear(P.env.epoch_utc));
    end
else
    % Reproducible default
    dy = 2026.0;
end

gen = 14;
if isfield(P, 'env') && isfield(P.env, 'igrf_generation') && ~isempty(P.env.igrf_generation)
    gen = double(P.env.igrf_generation);
end

XYZ_nT = igrfmagm(height_km, lat_deg, lon_deg, dy, gen); % [north east down] in nT
Bn = XYZ_nT(1);
Be = XYZ_nT(2);
Bd = XYZ_nT(3);

% NED -> ECEF (still nT)
[Bx_ecef, By_ecef, Bz_ecef] = ned2ecefv(Bn, Be, Bd, lat_deg, lon_deg, 'degrees');
B_ecef_nT = [Bx_ecef; By_ecef; Bz_ecef];

% ECEF -> ECI
B_eci_nT = R3.' * B_ecef_nT;

% nT -> Tesla
B_eci_T = 1e-9 * B_eci_nT;
end

