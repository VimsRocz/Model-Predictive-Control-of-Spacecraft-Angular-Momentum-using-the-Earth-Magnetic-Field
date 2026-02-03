function B_body = orbit_propagator_simple(k, P)
%ORBIT_PROPAGATOR_SIMPLE Synthetic B-field vector in body frame.
% This proxy is enough to validate MPC logic without full orbit/IGRF.

t = (k-1)*P.Ts;

% Rotate a nominal field direction over time (mimics orbit evolution)
w_orb = 2*pi/(90*60); % ~90 min orbit rad/s
theta = w_orb * t;

% Simple varying inertial-ish field
B_inertial = P.B0_T * [cos(theta); 0.3*sin(theta); 0.6*cos(0.5*theta)];

% Assume body frame ~ inertial frame for momentum-only model
B_body = B_inertial;
end
