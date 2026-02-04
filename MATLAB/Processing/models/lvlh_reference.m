function [q_ib_des, w_des_body, C_bi_des] = lvlh_reference(r_eci_m, v_eci_ms)
%LVLH_REFERENCE LVLH (nadir) attitude reference from orbit r,v in ECI.
%
% Outputs
%   q_ib_des  : 1x4 quaternion (inertial->body), scalar-first
%   w_des_body: 3x1 desired angular velocity expressed in body [rad/s]
%   C_bi_des  : 3x3 DCM inertial->body
%
% Convention
% - Body axes are aligned with LVLH:
%     +X along (approximately) velocity direction
%     +Z points to nadir (toward Earth center)
%     +Y completes right-hand rule

r = r_eci_m(:);
v = v_eci_ms(:);

rn = norm(r);
vn = norm(v);
if rn < 1e-9 || vn < 1e-9
    error("lvlh_reference requires non-zero r and v.");
end

rhat = r / rn;
vhat = v / vn;

z_b_in_eci = -rhat;

% Make x orthogonal to z (robust for non-perfectly circular orbits)
x_tmp = vhat - dot(vhat, z_b_in_eci) * z_b_in_eci;
xn = norm(x_tmp);
if xn < 1e-12
    % Degenerate; pick an arbitrary orthogonal axis.
    a = [1;0;0];
    if abs(dot(a, z_b_in_eci)) > 0.9
        a = [0;1;0];
    end
    x_tmp = cross(a, z_b_in_eci);
    xn = norm(x_tmp);
end
x_b_in_eci = x_tmp / xn;

y_b_in_eci = cross(z_b_in_eci, x_b_in_eci);
y_b_in_eci = y_b_in_eci / max(norm(y_b_in_eci), 1e-12);

% DCM inertial->body has rows equal to body axes expressed in inertial.
C_bi_des = [x_b_in_eci.'; y_b_in_eci.'; z_b_in_eci.'];

q_ib_des = dcm2quat(C_bi_des); % 1x4, scalar-first

% Orbital angular velocity in ECI (instantaneous)
w_orb_eci = cross(r, v) / max(rn^2, 1e-12);

% Desired angular velocity in body coordinates
w_des_body = C_bi_des * w_orb_eci;
end

