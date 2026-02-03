function tau_ext = ext_torques(k, P)
%EXT_TORQUES Simple external torque model (constant + small periodic)

t = (k-1)*P.Ts;
tau_ext = P.tau_ext_const + 1e-6*[sin(2*pi*t/600); cos(2*pi*t/800); sin(2*pi*t/1000)];
end
