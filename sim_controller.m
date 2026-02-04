function m = sim_controller(x, B, t)
%SIM_CONTROLLER Wrapper for Simulink Interpreted MATLAB Function block.

persistent P m_prev last_t
if isempty(P)
    P = evalin('base','P');
end

if isempty(last_t) || t < last_t
    m_prev = zeros(3,1);
end
last_t = t;

if isempty(m_prev)
    m_prev = zeros(3,1);
end

step = floor(t/P.Ts) + 1;

% Choose controller (default to baseline for reproducible MATLAB vs Simulink comparison)
mode = "baseline";
if isfield(P,'simulink') && isfield(P.simulink,'controller')
    mode = string(P.simulink.controller);
end
mode = lower(mode);

if mode == "baseline"
    m = baseline_mtq(x, B, P);
else
    Bprov   = @(kk) orbit_propagator_simple(kk, P);
    tauprov = @(kk) ext_torques(kk, P);
    m = mpc_mtq_qp(x, step, P, Bprov, tauprov, m_prev);
end

m_prev = m;
end
