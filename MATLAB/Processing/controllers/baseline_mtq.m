function m_cmd = baseline_mtq(x, B, P)
%BASELINE_MTQ Instantaneous baseline: tau_des = -K*x, m = (B x tau_des)/||B||^2

tau_des = -P.KH * x(:);

Bn2 = max(dot(B,B), 1e-12);
m_cmd = cross(B, tau_des) / Bn2;

% Saturate per axis
m_cmd = min(max(m_cmd, -P.m_max), P.m_max);
end
