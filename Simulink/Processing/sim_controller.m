function m = sim_controller(x, B, t)
%SIM_CONTROLLER Wrapper for Simulink Interpreted MATLAB Function block.

persistent P m_prev last_t last_run_id
[baseP, base_id] = load_base_P();
if isempty(P) || isempty(last_run_id) || ~strcmp(last_run_id, base_id)
    P = baseP;
    m_prev = [];
    last_t = [];
    last_run_id = base_id;
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

function [Pbase, run_id] = load_base_P()
Pbase = evalin('base','P');
run_id = "";
if isfield(Pbase,'simulink') && isfield(Pbase.simulink,'run_id')
    run_id = string(Pbase.simulink.run_id);
end
run_id = char(run_id);
end
