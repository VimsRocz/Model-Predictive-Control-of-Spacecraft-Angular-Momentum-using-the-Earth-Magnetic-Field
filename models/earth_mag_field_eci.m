function B_eci_T = earth_mag_field_eci(r_eci_m, P, t_s)
%EARTH_MAG_FIELD_ECI Earth magnetic field vector in ECI frame [Tesla].
%
% Select model via:
%   P.env.bfield_model = "igrf" | "dipole"
%
% Notes
% - "dipole" uses a simple tilted dipole model (fast, analytic).
% - "igrf" uses Aerospace Toolbox `igrfmagm` (IGRF-14) with a simple
%   Earth rotation ECI<->ECEF model (no precession/nutation).
%
% The returned vector is expressed in the inertial/ECI frame.

if nargin < 3
    t_s = 0;
end

model = "dipole";
if isfield(P, "env") && isfield(P.env, "bfield_model")
    model = string(P.env.bfield_model);
end
model = lower(model);

switch model
    case "dipole"
        B_eci_T = earth_dipole_field(r_eci_m, P, t_s);
    case "igrf"
        B_eci_T = earth_igrf_field(r_eci_m, P, t_s);
    otherwise
        error("Unknown bfield_model '%s'. Use 'dipole' or 'igrf'.", model);
end
end

